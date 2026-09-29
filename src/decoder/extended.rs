//! 12-bit DCT-based decoding (extended sequential and progressive processes).
//!
//! Only single-component images are supported. All quantized coefficients are
//! kept until the end of the image, then dequantized and transformed here with a
//! separable floating point IDCT that has a fixed evaluation order.
use crate::error::{Error, Result};
use crate::parser::{CodingProcess, FrameInfo};
use alloc::borrow::ToOwned;
use alloc::sync::Arc;
use alloc::vec;
use alloc::vec::Vec;

/// Returns whether `frame` uses a DCT-based process with 12-bit samples.
pub(crate) fn is_dct_12_bit(frame: &FrameInfo) -> bool {
    frame.precision == 12 && frame.coding_process != CodingProcess::Lossless
}

/// Returns whether all quantized coefficients of `frame` are kept until the end of the image.
pub(crate) fn stores_coefficients(frame: &FrameInfo) -> bool {
    frame.coding_process == CodingProcess::DctProgressive || is_dct_12_bit(frame)
}

/// Dequantizes and transforms the stored coefficients of a 12-bit single-component image.
///
/// Samples are returned as native-endian `u16` values, like lossless images with more than
/// 8 bits of precision.
pub(crate) fn compute_image_extended(
    frame: &FrameInfo,
    coefficients: &[Vec<i16>],
    quantization_tables: &[Option<Arc<[u16; 64]>>; 4],
) -> Result<Vec<u8>> {
    let component = match frame.components.as_slice() {
        [component] => component,
        _ => {
            return Err(Error::Format(
                "12-bit DCT decoding requires exactly one component".to_owned(),
            ))
        }
    };
    let quantization_table = quantization_tables[component.quantization_table_index]
        .as_ref()
        .ok_or_else(|| Error::Format("use of unset quantization table".to_owned()))?;

    let blocks_wide = usize::from(component.block_size.width);
    let blocks_high = usize::from(component.block_size.height);
    let coefficients = match coefficients.first() {
        Some(coefficients) if coefficients.len() == blocks_wide * blocks_high * 64 => coefficients,
        _ => return Err(Error::Format("not all components have data".to_owned())),
    };

    let width = usize::from(frame.output_size.width);
    let height = usize::from(frame.output_size.height);
    if width > blocks_wide * 8 || height > blocks_high * 8 {
        return Err(Error::Format(
            "12-bit image is larger than its coded blocks".to_owned(),
        ));
    }

    let basis = idct_basis();
    let level_shift = f64::from(1u16 << (frame.precision - 1));
    let max_value = f64::from((1u16 << frame.precision) - 1);
    let mut samples = vec![0u16; width * height];

    for block_y in 0..blocks_high {
        let top = block_y * 8;
        if top >= height {
            break;
        }
        for block_x in 0..blocks_wide {
            let left = block_x * 8;
            if left >= width {
                break;
            }
            let offset = (block_y * blocks_wide + block_x) * 64;
            let block = idct_block(
                &coefficients[offset..offset + 64],
                quantization_table,
                &basis,
            );
            for y in 0..8.min(height - top) {
                for x in 0..8.min(width - left) {
                    let value = (block[y * 8 + x] + level_shift).round();
                    samples[(top + y) * width + left + x] = value.clamp(0.0, max_value) as u16;
                }
            }
        }
    }

    Ok(samples
        .iter()
        .flat_map(|sample| sample.to_ne_bytes())
        .collect())
}

/// `basis[x][u] = C(u) * cos((2x + 1) * u * pi / 16) / 2`, with `C(0) = 1 / sqrt(2)`.
fn idct_basis() -> [[f64; 8]; 8] {
    let mut basis = [[0.0; 8]; 8];
    for (x, row) in basis.iter_mut().enumerate() {
        for (u, value) in row.iter_mut().enumerate() {
            let scale = if u == 0 {
                core::f64::consts::FRAC_1_SQRT_2
            } else {
                1.0
            };
            let angle = ((2 * x + 1) * u) as f64 * core::f64::consts::PI / 16.0;
            *value = scale * angle.cos() / 2.0;
        }
    }
    basis
}

/// Inverse DCT of one block in natural (row-major) coefficient order, without level shift.
fn idct_block(
    coefficients: &[i16],
    quantization_table: &[u16; 64],
    basis: &[[f64; 8]; 8],
) -> [f64; 64] {
    let mut dequantized = [0.0f64; 64];
    for (value, (&coefficient, &quantizer)) in dequantized
        .iter_mut()
        .zip(coefficients.iter().zip(quantization_table.iter()))
    {
        *value = f64::from(i32::from(coefficient) * i32::from(quantizer));
    }

    // Rows: horizontal frequencies u to horizontal positions x.
    let mut rows = [0.0f64; 64];
    for v in 0..8 {
        for x in 0..8 {
            let mut sum = 0.0;
            for u in 0..8 {
                sum += dequantized[v * 8 + u] * basis[x][u];
            }
            rows[v * 8 + x] = sum;
        }
    }

    // Columns: vertical frequencies v to vertical positions y.
    let mut output = [0.0f64; 64];
    for y in 0..8 {
        for x in 0..8 {
            let mut sum = 0.0;
            for v in 0..8 {
                sum += rows[v * 8 + x] * basis[y][v];
            }
            output[y * 8 + x] = sum;
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dc_only_block_is_flat() {
        let mut coefficients = [0i16; 64];
        // A DC value of 8 * 100 reconstructs a flat block of 100 before level shift.
        coefficients[0] = 800;
        let block = idct_block(&coefficients, &[1; 64], &idct_basis());
        for value in block {
            assert!((value - 100.0).abs() < 1e-9, "{value}");
        }
    }

    #[test]
    fn single_ac_coefficient_matches_definition() {
        let mut coefficients = [0i16; 64];
        coefficients[1] = 64; // u = 1, v = 0
        let block = idct_block(&coefficients, &[2; 64], &idct_basis());
        for y in 0..8 {
            for x in 0..8 {
                let expected = 0.25
                    * core::f64::consts::FRAC_1_SQRT_2
                    * 128.0
                    * (((2 * x + 1) as f64) * core::f64::consts::PI / 16.0).cos();
                assert!((block[y * 8 + x] - expected).abs() < 1e-9);
            }
        }
    }
}

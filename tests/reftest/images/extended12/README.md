Synthetic 12-bit DCT fixtures, single component.

The `.jpg` files were encoded by libjpeg-turbo 3.1.4 `cjpeg`. The `.png` references were
decoded by its `djpeg -dct int -pnm` and stored as 16-bit grayscale with values 0-4095.
`reftest_extended_12_bit` allows a difference of one level for IDCT rounding.

Source `source.pgm` (37 x 23, maxval 4095):

`v(x, y) = (x*97 + y*53 + ((x*y) % 37)*29) % 4096`, then `v //= 3` where `(x//5 + y//4) % 2 == 0`.

| File | Command |
| --- | --- |
| `sequential_37x23_q90.jpg` | `cjpeg -precision 12 -quality 90 source.pgm` |
| `sequential_37x23_restart1.jpg` | `cjpeg -precision 12 -quality 90 -restart 1 source.pgm` |
| `sequential_37x23_q5_16bit_dqt.jpg` | `cjpeg -precision 12 -quality 5 source.pgm` (16-bit DQT) |
| `sequential_37x23_sampling2x2.jpg` | `cjpeg -precision 12 -quality 90 -sample 2x2 source.pgm` |
| `progressive_37x23_q90.jpg` | `cjpeg -precision 12 -quality 90 -progressive source.pgm` |
| `sequential_64x16_q100_checker.jpg` | `cjpeg -precision 12 -quality 100 checker.pgm` |

`checker.pgm` (64 x 16) alternates flat 8 x 8 blocks of 0 and 4095. That forces DC
difference categories up to 15.

`../../unsupported/extended12_rgb_37x23.jpg` is a 3-component 12-bit image from the same
source, which must return `UnsupportedFeature::SamplePrecision(12)`.

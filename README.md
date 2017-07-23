# Recbf [![Build Status](https://travis-ci.org/mthh/recbf.svg?branch=master)](https://travis-ci.org/mthh/recbf)

__Rust library for recursive bilateral filtering.__

Example usage with `rust-image` library for reading/writing image files :
```rust
let img = image::open(&Path::new("in.jpg")).unwrap();
let mut content = img.raw_pixels();
recursive_bf(&mut content, 0.03, 0.1, width, height, channel);
image::save_buffer(Path::new("out.jpg"), &content, width, height, img.color());
```

<table>
<tr>
<td><img src="https://raw.githubusercontent.com/mthh/recbf/master/examples/in.jpg" width="270px"><br/><p align="center">Original Image</p></td>
<td><img src="https://raw.githubusercontent.com/mthh/recbf/master/examples/out_verif.jpg" width="270"><br/><p align="center">Recursive Bilateral Filtering</p></td>
</tr>
</table>

Translated from [RecursiveBF](https://github.com/ufoym/RecursiveBF) C++ library.
For more details on the algorithm, please refer to the original paper :

```
@inproceedings{yang2012recursive,
    title={Recursive bilateral filtering},
    author={Yang, Qingxiong},
    booktitle={European Conference on Computer Vision},
    pages={399--413},
    year={2012},
    organization={Springer}
}
```

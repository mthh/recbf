extern crate libc;
use libc::{memcpy, c_void, size_t};
use std::mem::size_of;


pub fn recursive_bf(img : &mut [u8], sigma_spatial : f32, sigma_range : f32, width : u32, height : u32, channel : u32) {
    unsafe {  _recursive_bf_unsafe(img.as_mut_ptr(), sigma_spatial, sigma_range, width, height, channel) };
}

#[no_mangle]
pub unsafe fn _recursive_bf_unsafe(img : *mut u8, sigma_spatial : f32, sigma_range : f32, width : u32, height : u32, channel : u32) {
    let width_height : u32 = width * height;
    let width_channel : u32 = width * channel;
    let width_height_channel : u32 = width * height * channel;

    let mut b = vec![0.0f32; ((width_height_channel + width_height + width_channel + width) * 2) as usize];
    let  img_out_f : *mut f32 = b.as_mut_ptr();;
    let img_temp : *mut f32 = &mut *img_out_f.offset(width_height_channel as (isize)) as (*mut f32);
    let map_factor_a : *mut f32 = &mut *img_temp.offset(width_height_channel as (isize)) as (*mut f32);
    let map_factor_b : *mut f32 = &mut *map_factor_a.offset(width_height as (isize)) as (*mut f32);
    let slice_factor_a : *mut f32 = &mut *map_factor_b.offset(width_height as (isize)) as (*mut f32);
    let slice_factor_b : *mut f32 = &mut *slice_factor_a.offset(width_channel as (isize)) as (*mut f32);
    let line_factor_a : *mut f32 = &mut *slice_factor_b.offset(width_channel as (isize)) as (*mut f32);
    let line_factor_b : *mut f32 = &mut *line_factor_a.offset(width as (isize)) as (*mut f32);

    let mut range_table = Vec::with_capacity(256);
    let inv_sigma_range : f32 = 1.0f32 / (sigma_range * 255i32 as (f32));
    for i in 0..256i32 {
        range_table.push(((-i as (f32) * inv_sigma_range) as (f64)).exp() as (f32));
    }

    let mut alpha : f32 = (-(2.0f64).sqrt() / (sigma_spatial * width as (f32)) as (f64)).exp() as (f32);
    let mut ypr : f32;
    let mut ypg : f32;
    let mut ypb : f32;
    let mut ycr : f32;
    let mut ycg : f32;
    let mut ycb : f32;
    let mut fp : f32;
    let mut fc : f32;
    let mut inv_alpha_ : f32 = 1i32 as (f32) - alpha;
    for y in 0..height {
        let mut temp_x : *mut f32 = &mut *img_temp.offset((y * width_channel) as (isize)) as (*mut f32);
        let mut in_x : *mut u8 = &mut *img.offset((y * width_channel) as (isize)) as (*mut u8);
        let mut texture_x : *mut u8 = &mut *img.offset((y * width_channel) as (isize)) as (*mut u8);
        ypr = *in_x as f32;
        *temp_x = ypr;
        temp_x = temp_x.offset(1isize);
        in_x = in_x.offset(1isize);
        ypg = *in_x as f32;
        *temp_x = ypg;
        temp_x = temp_x.offset(1isize);
        in_x = in_x.offset(1isize);
        ypb = *in_x as f32;
        *temp_x = ypb;
        temp_x = temp_x.offset(1isize);
        in_x = in_x.offset(1isize);
        let mut tpr : u8 = *texture_x;
        texture_x = texture_x.offset(1isize);
        let mut tpg : u8 = *texture_x;
        texture_x = texture_x.offset(1isize);
        let mut tpb : u8 = *texture_x;
        texture_x = texture_x.offset(1isize);
        let mut temp_factor_x : *mut f32 = &mut *map_factor_a.offset((y * width) as (isize)) as (*mut f32);
        fp = 1f32;
        *temp_factor_x = fp;
        temp_factor_x = temp_factor_x.offset(1isize);
        for _ in 1..width {
            let tcr : u8 = *texture_x;
            texture_x = texture_x.offset(1isize);
            let tcg : u8 = *texture_x;
            texture_x = texture_x.offset(1isize);
            let tcb : u8 = *texture_x;
            texture_x = texture_x.offset(1isize);
            let dr : u8 = (tcr as (i32) - tpr as (i32)).abs() as (u8);
            let dg : u8 = (tcg as (i32) - tpg as (i32)).abs() as (u8);
            let db : u8 = (tcb as (i32) - tpb as (i32)).abs() as (u8);
            let range_dist : i32 = (dr as (i32) << 1i32) + dg as (i32) + db as (i32) >> 2i32;
            let weight : f32 = range_table[range_dist as (usize)];
            let alpha_ : f32 = weight * alpha;

            ycr = inv_alpha_ * (*in_x) as (f32) + alpha_ * ypr;
            in_x = in_x.offset(1isize);
            *temp_x = ycr;
            temp_x = temp_x.offset(1isize);

            ycg = inv_alpha_ * (*in_x) as (f32) + alpha_ * ypg;
            in_x = in_x.offset(1isize);
            *temp_x = ycg;
            temp_x = temp_x.offset(1isize);

            ycb = inv_alpha_ * (*in_x) as (f32) + alpha_ * ypb;
            in_x = in_x.offset(1isize);
            *temp_x = ycb;
            temp_x = temp_x.offset(1isize);

            tpr = tcr;
            tpg = tcg;
            tpb = tcb;
            ypr = ycr;
            ypg = ycg;
            ypb = ycb;
            fc = inv_alpha_ + alpha_ * fp;
            *temp_factor_x = fc;
            temp_factor_x = temp_factor_x.offset(1isize);
            fp = fc;
        }
        temp_x = temp_x.offset(-1isize);
        *temp_x = 0.5f32 * (*temp_x + *{ in_x = in_x.offset(-1isize); in_x } as (f32));
        temp_x = temp_x.offset(-1isize);
        *temp_x = 0.5f32 * (*temp_x + *{ in_x = in_x.offset(-1isize); in_x } as (f32));
        temp_x = temp_x.offset(-1isize);
        *temp_x = 0.5f32 * (*temp_x + *{ in_x = in_x.offset(-1isize);  in_x } as (f32));
        tpr = *{ texture_x = texture_x.offset(-1isize); texture_x };
        tpg = *{ texture_x = texture_x.offset(-1isize); texture_x };
        tpb = *{ texture_x = texture_x.offset(-1isize); texture_x };
        ypr = *in_x as (f32);
        ypg = *in_x as (f32);
        ypb = *in_x as (f32);
        temp_factor_x = temp_factor_x.offset(-1isize);
        *temp_factor_x = 0.5f32 * (*temp_factor_x + 1i32 as (f32));
        fp = 1f32;
        for _ in (0..width - 1u32).rev() {
            texture_x = texture_x.offset(-1isize);
            let tcr : u8 = *texture_x;
            texture_x = texture_x.offset(-1isize);
            let tcg : u8 = *texture_x;
            texture_x = texture_x.offset(-1isize);
            let tcb : u8 = *texture_x;
            let dr : u8 = (tcr as (i32) - tpr as (i32)).abs() as (u8);
            let dg : u8 = (tcg as (i32) - tpg as (i32)).abs() as (u8);
            let db : u8 = (tcb as (i32) - tpb as (i32)).abs() as (u8);
            let range_dist : i32 = (dr as (i32) << 1i32) + dg as (i32) + db as (i32) >> 2i32;
            let weight : f32 = range_table[range_dist as (usize)];
            let alpha_ : f32 = weight * alpha;
            in_x = in_x.offset(-1isize);
            ycr = inv_alpha_ * *in_x as (f32) + alpha_ * ypr;
            in_x = in_x.offset(-1isize);
            ycg = inv_alpha_ * *in_x as (f32) + alpha_ * ypg;
            in_x = in_x.offset(-1isize);
            ycb = inv_alpha_ * *in_x as (f32) + alpha_ * ypb;
            temp_x = temp_x.offset(-1isize);
            *temp_x = 0.5f32 * (*temp_x + ycr);
            temp_x = temp_x.offset(-1isize);
            *temp_x = 0.5f32 * (*temp_x + ycg);
            temp_x = temp_x.offset(-1isize);
            *temp_x = 0.5f32 * (*temp_x + ycb);
            tpr = tcr;
            tpg = tcg;
            tpb = tcb;
            ypr = ycr;
            ypg = ycg;
            ypb = ycb;
            fc = inv_alpha_ + alpha_ * fp;
            temp_factor_x = temp_factor_x.offset(-1isize);
            *temp_factor_x = 0.5f32 * (*temp_factor_x + fc);
            fp = fc;
        }
    }
    alpha = (-(2.0f64).sqrt() / (sigma_spatial * height as (f32)) as (f64)).exp() as (f32);
    inv_alpha_ = 1i32 as (f32) - alpha;
    let mut ycy : *mut f32;
    let mut ypy : *mut f32;
    let mut xcy : *mut f32;
    let mut tcy : *mut u8;
    let mut tpy : *mut u8;
    memcpy(
        img_out_f as (*mut c_void),
        img_temp as (*const c_void),
        size_of::<f32>().wrapping_mul(width_channel as usize) as size_t
    );
    let in_factor : *mut f32 = map_factor_a;
    let mut ycf : *mut f32;
    let mut ypf : *mut f32;
    let mut xcf : *mut f32;
    memcpy(
        map_factor_b as (*mut c_void), in_factor as (*const c_void),
        size_of::<f32>().wrapping_mul(width as usize) as size_t
    );
    for y in 1..height {
        tpy = &mut *img.offset(((y - 1u32) * width_channel) as (isize)) as (*mut u8);
        tcy = &mut *img.offset((y * width_channel) as (isize)) as (*mut u8);
        xcy = &mut *img_temp.offset((y * width_channel) as (isize)) as (*mut f32);
        ypy = &mut *img_out_f.offset(((y - 1u32) * width_channel) as (isize)) as (*mut f32);
        ycy = &mut *img_out_f.offset((y * width_channel) as (isize)) as (*mut f32);
        xcf = &mut *in_factor.offset((y * width) as (isize)) as (*mut f32);
        ypf = &mut *map_factor_b.offset(((y - 1u32) * width) as (isize)) as (*mut f32);
        ycf = &mut *map_factor_b.offset((y * width) as (isize)) as (*mut f32);
        for _ in 0..width {
            let dr : u8 = (*tcy as i32 - *tpy as i32).abs() as u8;
            tcy = tcy.offset(1isize);
            tpy = tpy.offset(1isize);
            let dg : u8 = (*tcy as i32 - *tpy as i32).abs() as u8;
            tcy = tcy.offset(1isize);
            tpy = tpy.offset(1isize);
            let db : u8 = (*tcy as i32 - *tpy as i32).abs() as u8;
            tcy = tcy.offset(1isize);
            tpy = tpy.offset(1isize);
            let range_dist : i32 = (dr as (i32) << 1i32) + dg as (i32) + db as (i32) >> 2i32;
            let weight : f32 = range_table[range_dist as (usize)];
            let alpha_ : f32 = weight * alpha;
            for _ in 0..channel {
                *ycy = inv_alpha_ * (*xcy) + alpha_ * (*ypy);
                ypy = ypy.offset(1isize);
                xcy = xcy.offset(1isize);
                ycy = ycy.offset(1isize);
            }
            *ycf = inv_alpha_ * (*xcf) + alpha_ * (*ypf);
            ycf = ycf.offset(1isize);
            xcf = xcf.offset(1isize);
            ypf = ypf.offset(1isize);
        }
    }
    let h1 : u32 = height - 1u32;
    ycf = line_factor_a;
    ypf = line_factor_b;
    memcpy(
        ypf as (*mut c_void),
        &mut *in_factor.offset( (h1 * width) as (isize)) as (*mut f32) as (*const c_void),
        size_of::<f32>().wrapping_mul(width as usize) as size_t
    );
    for x in 0..width {
        *map_factor_b.offset(
             (h1 * width + x) as (isize)
         ) = 0.5f32 * (*map_factor_b.offset((h1 * width + x) as (isize)) + *ypf.offset(x as (isize)));
    }
    ycy = slice_factor_a;
    ypy = slice_factor_b;
    memcpy(
        ypy as (*mut c_void),
        &mut *img_temp.offset( (h1 * width_channel) as (isize)) as (*mut f32) as (*const c_void),
        size_of::<f32>().wrapping_mul(width_channel as usize) as size_t
    );
    let mut k : i32 = 0i32;
    for x in 0..width {
        for c in 0..channel {
            let idx : isize = ((h1 * width + x) * channel + c) as isize;
            *img_out_f.offset(idx) = 0.5f32 * (
                    *img_out_f.offset(idx) + *ypy.offset(k as (isize))
                ) / *map_factor_b.offset((h1 * width + x) as (isize));
            k += 1;
        }
    }
    for y in (0..h1).rev() {
        tpy = &mut *img.offset( ((y + 1u32) * width_channel) as (isize) ) as (*mut u8);
        tcy = &mut *img.offset( (y * width_channel) as (isize) ) as (*mut u8);
        xcy = &mut *img_temp.offset( (y * width_channel) as (isize) ) as (*mut f32);
        xcf = &mut *in_factor.offset((y * width) as (isize)) as (*mut f32);
        let mut ycy_ : *mut f32 = ycy;
        let mut ypy_ : *mut f32 = ypy;
        let mut out_ : *mut f32 = &mut *img_out_f.offset( (y * width_channel) as (isize) ) as (*mut f32);
        let mut ycf_ : *mut f32 = ycf;
        let mut ypf_ : *mut f32 = ypf;
        let mut factor_ : *mut f32 = &mut *map_factor_b.offset((y * width) as (isize)) as (*mut f32);
        for _ in 0..width {
            let dr = (*tcy as i32 - *tpy as i32).abs() as (u8);
            tcy = tcy.offset(1isize);
            tpy = tpy.offset(1isize);
            let dg = (*tcy as i32 - *tpy as i32).abs() as (u8);
            tcy = tcy.offset(1isize);
            tpy = tpy.offset(1isize);
            let db = (*tcy as i32 - *tpy as i32).abs() as (u8);
            tcy = tcy.offset(1isize);
            tpy = tpy.offset(1isize);

            let range_dist : i32 = (dr as (i32) << 1i32) + dg as (i32) + db as (i32) >> 2i32;
            let weight : f32 = range_table[range_dist as (usize)];
            let alpha_ : f32 = weight * alpha;
            let fcc = inv_alpha_ * (*xcf) + alpha_ * (*ypf_);
            xcf = xcf.offset(1isize);
            ypf_ = ypf_.offset(1isize);
            *ycf_ = fcc;
            ycf_ = ycf_.offset(1isize);
            *factor_ = 0.5f32 * (*factor_ + fcc);
            for _ in 0..channel {
                let ycc = inv_alpha_ * (*xcy) + alpha_ * (*ypy_);
                xcy = xcy.offset(1isize);
                ypy_ = ypy_.offset(1isize);
                *ycy_ = ycc;
                ycy_ = ycy_.offset(1isize);
                *out_ = 0.5f32 * ((*out_) + ycc) / (*factor_);
                out_ = out_.offset(1isize);
            }
            factor_ = factor_.offset(1isize);
        }
        memcpy(
            ypy as (*mut c_void), ycy as (*const c_void),
            size_of::<f32>().wrapping_mul(width_channel as usize) as size_t
        );
        memcpy(
            ypf as (*mut c_void), ycf as (*const c_void),
            size_of::<f32>().wrapping_mul(width as usize) as size_t
        );
    }
    for i in 0..width_height_channel as isize {
        *img.offset(i) = *img_out_f.offset(i) as (u8);
    }
}

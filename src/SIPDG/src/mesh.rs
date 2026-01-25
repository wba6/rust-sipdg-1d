use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int};
use std::ptr;
use gmsh_sys::*;

/// Generates a 1D mesh using raw gmsh-sys C-bindings.
pub fn generate_1d_gmesh(start: f64, end: f64, num_elements: usize) -> Vec<f64> {
    unsafe {
        // Initialize with null pointers for error handling
        gmshInitialize(0, ptr::null_mut(), 1, ptr::null_mut());

        let model_name = CString::new("SIPDG_Model").unwrap();
        gmshModelAdd(model_name.as_ptr() as *const c_char, ptr::null_mut());

        let lc = (end - start) / (num_elements as f64);
        
        // Define geometry
        gmshModelGeoAddPoint(start as c_double, 0.0, 0.0, lc as c_double, 1, ptr::null_mut());
        gmshModelGeoAddPoint(end as c_double, 0.0, 0.0, lc as c_double, 2, ptr::null_mut());
        gmshModelGeoAddLine(1, 2, 1, ptr::null_mut());
        gmshModelGeoSynchronize(ptr::null_mut());

        // Generate 1D Mesh
        gmshModelMeshGenerate(1, ptr::null_mut());

        // Prepare pointers for out-parameters
        let mut node_tags_ptr: *mut usize = ptr::null_mut();
        let mut node_tags_len: usize = 0;
        let mut coords_ptr: *mut c_double = ptr::null_mut();
        let mut coords_len: usize = 0;
        let mut params_ptr: *mut c_double = ptr::null_mut();
        let mut params_len: usize = 0;

        // FIXED: Reordered arguments to match gmsh-sys 0.1.2 bindings
        // The pattern is: (out_tags, out_coords, out_params, dim, tag, includeBoundary, returnParametric, ierr)
        gmshModelMeshGetNodes(
            &mut node_tags_ptr, 
            &mut node_tags_len,
            &mut coords_ptr, 
            &mut coords_len,
            &mut params_ptr, 
            &mut params_len,
            -1 as c_int,     // dim
            -1 as c_int,     // tag
            1 as c_int,      // includeBoundary (true)
            0 as c_int,      // returnParametric (false)
            ptr::null_mut()  // ierr
        );

        // Convert the raw C buffer into a Rust slice and collect X coordinates
        let raw_coords = std::slice::from_raw_parts(coords_ptr, coords_len);
        let mut x_coords: Vec<f64> = raw_coords
            .chunks(3)
            .map(|c| c[0])
            .collect();

        // Clean up Gmsh-allocated memory
        gmshFree(node_tags_ptr as *mut _);
        gmshFree(coords_ptr as *mut _);
        gmshFree(params_ptr as *mut _);

        gmshFinalize(ptr::null_mut());

        // Sort for SIPDG element ordering
        x_coords.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        x_coords
    }
}

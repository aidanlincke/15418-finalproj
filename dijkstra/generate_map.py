from pyproj import Transformer
from rasterio.windows import from_bounds
import rasterio
import numpy as np
import struct
from skimage.transform import resize
from affine import Affine

def world_to_pixel(transform, x, y):
    col = int((x - transform.c) / transform.a)
    row = int((y - transform.f) / transform.e)
    return row, col

if __name__ == "__main__":
    with rasterio.open("problems/dc.tif") as tif:
        transformer = Transformer.from_crs("EPSG:4326", tif.crs, always_xy=True)

        wm_x, wm_y = transformer.transform(-77.0353, 38.8895) # Washington Monument
        dca_x, dca_y = transformer.transform(-77.0402, 38.8512) # DCA

        x_min = min(wm_x, dca_x) - 100
        x_max = max(wm_x, dca_x) + 100
        y_min = min(wm_y, dca_y) - 100
        y_max = max(wm_y, dca_y) + 100

        window = from_bounds(x_min, y_min, x_max, y_max, tif.transform)
        cropped = tif.read(1, window=window)

    cropped = np.nan_to_num(cropped, nan=0.0).astype(np.float32)
    cropped -= np.min(cropped) # Normalize to >=0

    original_rows, original_cols = cropped.shape
    target_rows = 3000
    target_cols = int(original_cols * (target_rows / original_rows))

    resized = resize(
        cropped, 
        (target_rows, target_cols), 
        order=1, 
        preserve_range=True, 
        anti_aliasing=True
    ).astype(np.float32)

    cropped = resized
    rows, cols = cropped.shape

    scale_x = (x_max - x_min) / cols
    scale_y = (y_max - y_min) / rows
    transform = Affine.translation(x_min, y_min) * Affine.scale(scale_x, scale_y)

    start_row, start_col = world_to_pixel(transform, dca_x, dca_y)
    goal_row, goal_col = world_to_pixel(transform, wm_x, wm_y)

    with open("problems/dc.bin", "wb") as f:
        f.write(struct.pack("ii", rows, cols))
        f.write(struct.pack("ii", start_row, start_col))
        f.write(struct.pack("ii", goal_row, goal_col))
        f.write(cropped.tobytes())
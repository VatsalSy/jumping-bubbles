#import subprocess as sp
import numpy as np
import pyvista as pv
import os
import multiprocessing
import gc  # Garbage collection

def parse_vertex(vertex_str):
    return np.fromstring(vertex_str, sep=' ')

def create_polydata(points, lines_array, faces_array):
    polydata = pv.PolyData()
    polydata.points = np.array(points)
    polydata.lines = lines_array
    polydata.faces = faces_array
    return polydata

def run_process(command):
    p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    return stderr.decode("utf-8").strip()

def process_cells(cell_data):
    lines_array = np.array([[2, 0, 1], [2, 1, 2], [2, 2, 3], [2, 3, 0]], dtype=np.int32)
    faces_array = np.array([[4, 0, 1, 2, 3]], dtype=np.int32)
    cells = [[parse_vertex(point) for point in cell.split('\n')] for cell in cell_data.split('\n\n')]
    return [create_polydata(cell, lines_array, faces_array) for cell in cells]

def process_facets(facet_data):
    vertex_to_index = {}
    points, cells = [], []
    for facet in facet_data.split('\n\n'):
        cell = []
        for vertex_str in facet.split('\n'):
            vertex = parse_vertex(vertex_str)
            vertex_tuple = tuple(vertex)  # Convert numpy array to tuple
            if vertex_tuple not in vertex_to_index:
                vertex_to_index[vertex_tuple] = len(points)
                points.append(vertex)
            cell.append(vertex_to_index[vertex_tuple])
        cells.extend([len(cell)] + cell)
    return pv.PolyData(np.array(points, dtype=np.float64), faces=np.array(cells, dtype=np.int32))

def reflect_mesh(mesh, normals):
    for normal in normals.values():
        mesh = mesh.merge(mesh.reflect(normal))
    return mesh

def gettingGrid(filename):
    cell_data = run_process(["./getCells_bottomPlate", filename])
    cells = process_cells(cell_data)
    mesh = pv.MultiBlock(cells).combine()
    return reflect_mesh(mesh, {'x': [1, 0, 0], 'z': [0, 0, 1]})

def gettingFacets3D(filename):
    facet_data = run_process(["./getFacets3D", filename])
    mesh = process_facets(facet_data)
    return reflect_mesh(mesh, {'x': [1, 0, 0], 'z': [0, 0, 1]})

def process_and_save_image(t, base_filename, image_folder):
    
    filename = base_filename % t
    image_filename = f"{image_folder}/{int(t*1e4):06d}.png"
    
    if not os.path.exists(filename):
        print(f"File {filename} does not exist")
        return
    if os.path.exists(image_filename):
        print(f"Image {image_filename} already exists")
        return
    
    print(f"Processing {t}")

    cells = gettingGrid(filename)
    poly_data = gettingFacets3D(filename)

    plotter = pv.Plotter(off_screen=True)

    try: 
        plotter.add_mesh(cells, color='grey') 
        plotter.add_mesh(cells, color='black', style='wireframe')
        plotter.add_mesh(poly_data, color='orange')
        camera_position = (7.75, 3.50, -2.50)
        focal_point = (0.0, 0.0, 0.0)
        up_vector = (-0.30, 0.95, 0.025)
        CameraParallelScale = 2.2

        plotter.camera.position = camera_position
        plotter.camera.focal_point = focal_point
        plotter.camera.up = up_vector
        plotter.camera.parallel_scale = CameraParallelScale

        plotter.add_text(f"t = {t:.2f}", position='upper_right', font_size=15, color='black', font_file='cmunrm.ttf')
        plotter.screenshot(image_filename)
        plotter.close()
    except:
        print(f"Error in {t}")
    finally:                  
        # Explicitly deleting large objects and calling garbage collection
        del cells, poly_data, plotter
        gc.collect()

    # # Enable off-screen rendering
    # plotter = pv.Plotter(off_screen=True)
    # plotter.add_mesh(cells, color='grey') 
    # plotter.add_mesh(cells, color='black', style='wireframe')
    # plotter.add_mesh(poly_data, color='orange')

    # camera_position = (7.75, 3.50, -2.50)
    # focal_point = (0.0, 0.0, 0.0)
    # up_vector = (-0.30, 0.95, 0.025)
    # CameraParallelScale = 2.2

    # plotter.camera.position = camera_position
    # plotter.camera.focal_point = focal_point
    # plotter.camera.up = up_vector
    # plotter.camera.parallel_scale = CameraParallelScale

    # plotter.add_text(f"t = {t:.2f}", position='upper_right', font_size=15, color='black', font_file='cmunrm.ttf')

    # Render and save the screenshot
    # plotter.show(auto_close=False)
    # plotter.screenshot(image_filename)
    # plotter.close()

def process_single_time_step(args):
    process_and_save_image(*args)
    gc.collect()  # Additional garbage collection after each process

if __name__ == '__main__':
    base_filename = "intermediate/snapshot-%5.4f"
    image_folder = "Video"
    
    if not os.path.exists(image_folder):
        os.makedirs(image_folder)

    time_step = 0.01
    end_time = 8
    time_steps = np.arange(0, end_time + time_step, time_step)

    args = [(t, base_filename, image_folder) for t in time_steps]

    with multiprocessing.Pool(processes=2) as pool:
        pool.map(process_single_time_step, args)

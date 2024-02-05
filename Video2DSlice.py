import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import StrMethodFormatter
import concurrent.futures
import multiprocessing
import argparse

# Set matplotlib parameters
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
AxesLabel, TickLabel = [50, 20]

def get_data(exe):
    result = sp.run(exe, capture_output=True, text=True)
    data = np.fromstring(result.stderr, sep=' ').reshape((-1, 5))
    return data.T

def plot_subplot(ax, data, ymin, xmax, ymax):
    y, x, f, vel, D2 = data

    ax.tricontour(x, y, f, levels=[0.5], colors='green', linewidths=5)
    ax.tricontour(-x, y, f, levels=[0.5], colors='green', linewidths=5)
    cntrl1 = ax.tricontourf(-x, y, vel, levels=np.linspace(0, 5, 500), cmap='Purples', extend='max')
    cntrl2 = ax.tricontourf(x, y, D2, levels=np.linspace(-2, 2, 400), cmap='hot_r', extend='both')

    ax.plot([0, 0], [ymin, ymax], '--', color='grey', linewidth=4)

    rect = patches.Rectangle((-xmax, ymin), 2*xmax, ymax-ymin, linewidth=6, edgecolor='k', facecolor='none')
    ax.add_patch(rect)

    ax.set_aspect('equal')
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis('off')
    return cntrl1, cntrl2

def plot_data(filename, ImageName, ymin, xmax, ymax, n, Oh, LINEAR, t):
    # Adjust the figsize to accommodate both subplots with their aspect ratios
    fig_width = 24  # Adjust this as needed
    fig_height = fig_width * (ymax - ymin) / xmax  
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Create subplots
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # Execute commands to get data
    commands = [
        ["./getDataXSlice", filename, str(ymax), str(xmax), str(0.), str(n), str(Oh), LINEAR],
        ["./getDataZSlice", filename, str(ymax), str(xmax), str(0.), str(n), str(Oh), LINEAR]
    ]
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = executor.map(get_data, commands)

    # Plot subplots with adjusted aspect ratios
    for ax, data in zip([ax1, ax2], results):
        cntrl1, cntrl2 = plot_subplot(ax, data, ymin, xmax, ymax)
        ax.set_aspect('equal')  # Enforce equal aspect ratio

    # Manually adjust positions of the subplots
    # Note: Adjust these values as needed to align the subplots
    ax1.set_position([0.095, 0.1, 0.4, 0.8])
    ax2.set_position([0.4975, 0.1, 0.4, 0.8])

    # add colorbars for cntrl1 to the left of ax1 and cntrl2 to the right of ax2
    l, b, w, h = ax1.get_position().bounds
    cb1 = fig.add_axes([l+0.05*w, b-0.075*h, 0.9*w, 0.01])
    c1 = plt.colorbar(cntrl1,cax=cb1,orientation='horizontal')
    c1.set_label(r'$\|u_i\|/\sqrt{\gamma/\rho_lR_0}$',fontsize=TickLabel, labelpad=5)
    c1.ax.tick_params(labelsize=TickLabel)
    c1.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))

    l, b, w, h = ax2.get_position().bounds
    cb2 = fig.add_axes([l+0.05*w, b-0.075*h, 0.9*w, 0.01])
    c2 = plt.colorbar(cntrl2,cax=cb2,orientation='horizontal')
    c2.set_label(r'$\log_{10}\left(2Oh\mathcal{D}_{ij}\mathcal{D}_{ij}\right)$',fontsize=TickLabel)
    c2.ax.tick_params(labelsize=TickLabel)
    c2.ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    
    ax1.set_title(r'$t/\sqrt{\rho_lR_0^3/\gamma} = %3.2f$' % t, fontsize=TickLabel+10, pad=10)

    plt.savefig(ImageName, bbox_inches='tight')
    plt.close()
    # plt.show()

def process_single_time_step(t, base_filename, image_folder, ymin, Oh):
    filename = base_filename % t
    ImageName = os.path.join(image_folder, f'snapshot-{t:.4f}.png')
    
    if not os.path.exists(filename):
        print(f"File {filename} does not exist")
        return
    if os.path.exists(ImageName):
        print(f"Image {ImageName} already exists")
        return
    
    print(f"Processing {t}")

    xmax, ymax, n  = 2.5, 2.5, 256
    LINEAR = 'false'
    plot_data(filename, ImageName, ymin, xmax, ymax, n, Oh, LINEAR, t)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process fluid dynamics data.")
    parser.add_argument('--Oh', type=float, default=0.01, help='Ohnesorge number')

    args = parser.parse_args()

    base_filename = "intermediate/snapshot-%5.4f"
    image_folder = "Video2DSlice"
    ymin = -1.025
    Oh = args.Oh  # Retrieve Oh value from command line
    
    if not os.path.exists(image_folder):
        os.makedirs(image_folder)

    time_step = 0.01
    end_time = 8
    time_steps = np.arange(0, end_time + time_step, time_step)

    with multiprocessing.Pool(processes=4) as pool:
        pool.starmap(process_single_time_step, [(t, base_filename, image_folder, ymin, Oh) for t in time_steps])

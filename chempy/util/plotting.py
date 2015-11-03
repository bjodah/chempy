def mpl_outside_legend(ax, **kwargs):
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), **kwargs)

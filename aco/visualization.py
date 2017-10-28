import matplotlib.pyplot as plt

def prepare_route(positions, route):
    edges_x, edges_y = [],[]
    for position in route:
        x, y = positions[position]
        edges_x += [x]
        edges_y += [y]
    return edges_x, edges_y

def prepare_nodes( nodes ):
    return zip(*nodes)


class VTSP(object):
    """Visualize TSP"""
    def __init__(self, positions, invert):
        self.fig = plt.figure()
        self.ax1 = self.fig.add_subplot(111)
        self.invert = invert
        self.positions = positions

    def init(self):
        lon, lat = prepare_nodes(self.positions)
        if self.invert:
            lat, lon = lon, lat
        self.ax1.scatter(lon, lat, s=50, facecolors='none', edgecolors='r')
        self.route, = self.ax1.plot([],[], linestyle='-', color='y')
        plt.axis('off')
        plt.show(block=False)


    def update_route(self, route):

        route_x, route_y = prepare_route(self.positions, route)

        if self.invert:
            route_y, route_x = route_x, route_y

        self.route.set_xdata(route_x)
        self.route.set_ydata(route_y)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def freeze(self, best_distance):
        """Avoids closing the window"""
        plt.title('Best route: %.2f ' %(best_distance ))
        plt.show()
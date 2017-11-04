import numpy as np
import argparse
import math
from visualization import VTSP






def read_data( fname ):
    """Reads nodes from text file"""
    nodes = []
    positions = []
    with open(fname) as f:
        for line in f:
            parts = line.rstrip().split(' ')
            if len(parts) == 3 and parts[0].isdigit():
                num_node = int(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                # data starts with 1
                nodes.append(num_node - 1)
                positions.append((x,y))
    return nodes, positions


def euclidean(a, b):
    """Calculate the euclidean distance for two points"""
    return np.linalg.norm(a - b)


def compute_distances( nodes , positions):
    """Compute the euclidean distances among all the nodes"""
    distances = {}
    for a in nodes:
        for b in nodes:
            if a != b:
                # only stores a,b or b,a
                if not (a, b) in distances and not (b, a) in distances: # use a=distances((a,b)) or distances((b,a))
                    pair = (a, b)
                    distances[pair] = euclidean(positions[a], positions[b])
    return distances


class ACO(object):
    def __init__(self, nodes, positions):
        """Initialize the ACO"""
        self.nodes = np.array(nodes, dtype=int)
        self.positions = np.array(positions)
        self.distances = compute_distances(self.nodes, self.positions)
        self.pairs = { key: None for key in self.distances.keys() }

    def setup(self, alpha=None, beta=None, evaporation=None, num_ants=None, visualization=False):
        """Setup parameters and initializations"""
        self.alpha = alpha or 1.0
        self.beta = beta or 1.0
        self.evaporation = evaporation or 1.0
        self.num_ants = num_ants or 10
        # set all pheromones to 1.0
        self.pheromones = {key: 1.0 for key in self.distances.keys() }

        # this function is now optimized to be ran over vectors
        # vectorize __compute_probs__
        # Then we can use v_compute_probs with a numpy array
        self.v_compute_probs = np.vectorize(ACO.__compute_probs__)

        # initializes the visualization class
        self.visualization = visualization
        if self.visualization:
            self.vis = VTSP(self.positions, invert=True)

    def __correct_pair__(self, a, b):
        """Returns the correct pair associated to those nodes"""
        if (a, b) in self.pairs:
            return (a, b)
        return (b, a)

    @staticmethod
    ## Static method which is going to be optimized to run over a vector
    def __compute_probs__(x ,y, pheromones, distances, alpha, beta,  pick_pair):
        """Computes tau^alpha * (1/eta)^beta"""
        return math.pow(pheromones[pick_pair(x,y)], alpha) * math.pow(1.0 / distances[pick_pair(x,y)], beta)

    def run(self):
        """Runs everything"""
        if self.visualization:
           self.vis.init()

        best_distance = np.inf
        prev_best_distance = np.inf
        best_route = None
        min_convergence_iterations = self.num_ants
        convergence_iterations = min_convergence_iterations

        iter = 0
        while convergence_iterations:

            ## construct solutions
            routes, route_distances = self.__step__(self.num_ants, self.nodes, self.distances,
                                                    self.pheromones, self.alpha, self.beta)

            ## evaluate solutions
            i_min_distance = route_distances.min()
            if i_min_distance < best_distance:
                best_route = routes[np.argmin(route_distances)]
                best_distance = i_min_distance


            # if current best solution equals previous
            if prev_best_distance == best_distance:
                convergence_iterations -= 1
            else:
                convergence_iterations = min_convergence_iterations

            ## update pheromones
            self.__update_pheromones__(routes, route_distances, self.pheromones, self.evaporation)

            ## shows visualization each iterations
            ## otherwise the window becomes buggy
            if self.visualization:
                self.vis.update_route(best_route)

            print """Interation: %i Best solution distance: %f  %i""" % (iter, best_distance, convergence_iterations)
            iter += 1
            prev_best_distance = best_distance

        # after everything, this fixes the window
        if self.visualization:
            self.vis.freeze( best_distance )

        return best_distance, best_route

    def __step__(self, ants, nodes, distances, pheromones, alpha, beta):
        """Perform a complete step of the algorithm"""
        route_distances = np.zeros(ants)
        routes = np.zeros((ants, nodes.shape[0] + 1), dtype=int)
        # one computation per ant
        for i in range(self.num_ants):
            routes[i,:], route_distances[i] = self.__ant_solve__(nodes, distances, pheromones, alpha, beta)
        return routes, route_distances

    def __ant_solve__(self, nodes, distances, pheromones, alpha, beta):
        """Performs one complete solution"""
        remaining = nodes.shape[0] # num of remaining nodes
        rem_nodes = np.ones(nodes.shape, dtype=bool) # remaining nodes
        route = [] # route performed
        route_distance = 0.0
        init = np.random.choice(nodes) # initial position

        # init of the algorithm
        rem_nodes[init] = 0
        route.append(init)
        current = init
        remaining -= 1

        while remaining:
            # collects all remaining nodes
            rem_neighbors = nodes[rem_nodes]

            if len(rem_neighbors) > 1:
                next_node = self.__next_node__(current, rem_neighbors, distances, pheromones, alpha, beta )
            else:
                next_node = rem_neighbors[0]

            route.append(next_node)
            rem_nodes[next_node] = 0
            route_distance += distances[self.__correct_pair__(current, next_node)]
            remaining -= 1
            current = next_node

        # appends the initial node to complete the circular route and updates the distance
        route.append(init)
        route_distance += distances[self.__correct_pair__(current,  init)]

        return route, route_distance

    def __next_node__(self, _from, rem_neighbors, distances, pheromones, alpha, beta):
        """returns the next node based on the probabilities"""
        # applies the vectorized function to the rem_neighbors array
        rem_probs = self.v_compute_probs(rem_neighbors, _from,  pheromones, distances, alpha, beta,  self.__correct_pair__)
        sum_probs = rem_probs.sum() # sums everything to normalize it afterwards

        if sum_probs != 0.0:
            # picks one neighbor based on the probabilities normalized
            return rem_neighbors[ np.where(np.random.multinomial(1, rem_probs / sum_probs))][0]

        #if the sum of the probs is 0, just pick the argument with the maximum probability
        return rem_neighbors[ np.argmax(rem_probs)][0]

    def __update_pheromones__(self, routes, route_distances,  pheromones, evaporation):
        # apply evaporation to al edges
        for edge in pheromones.keys():
            pheromones[edge] *= (1. - evaporation)

        # update  pheromones based on the routes
        for i, route in enumerate(routes):
            for j in range(1, route.shape[0]):
                a = route[j]
                b = route[j - 1]
                # the longer the route, the lower the increase
                pheromones[self.__correct_pair__(a,b)] +=  1. / route_distances[i]



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, type=str, default=None,
                        help="TSP file containing the data")
    parser.add_argument("-a", "--alpha", required=False, type=float, default=1.0, help="Alpha value should be >= 0")
    parser.add_argument("-b", "--beta", required=False, type=float, default=2.0, help="Beta value should be >= 1")
    parser.add_argument("-e", "--evaporation", required=False, type=float, default=0.5, help="Evaporation value 0-1")
    parser.add_argument("-n", "--ants", required=False, type=int, default=15, help="Num of ants")
    parser.add_argument("-v", "--visualization", required=False, action='store_true', help="Visualize routes on real time")

    args = parser.parse_args()

    nodes, positions = read_data(args.file)


    aco = ACO(nodes, positions)
    aco.setup(args.alpha, args.beta, args.evaporation, args.ants, visualization=args.visualization)
    best_distance, best_route = aco.run()


    print """Overall best solution distance: %f
    Best route %s
    """ % (best_distance, " -> ".join([str(x + 1) for x in best_route]))
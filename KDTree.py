import numpy
class Node(object): 
    """
    A simply node for a binary k-d tree
    """
    def __init__(self,datum,left_child,right_child,axis=0,parent=None):
        """
        Create a new binary node
        """
        self.datum = datum
        self.axis = axis
        self.left = left_child
        self.right = right_child
        self.parent = parent

    @staticmethod 
    def inorder(node):
        """
        N.(Node) -> str

        Give a string representation of the in order traversal of a binary tree
        of nodes, with each node's contents separated by each others by spaces
        """
        if not(node.left or node.right):
            return str(node.datum)
        output = str(node.datum)
        if node.left:
            output = str(Node.inorder(node.left)) +" "+ output
        if node.right:
            output = output+" "+str(Node.inorder(node.right))
        return output

class KDTree(object):
    """
    A k-d tree equipped with a nearest neighbour method.
    """
    def __init__(self, points):
        """
        Create a new k-d tree with the given points
        """
        # TODO: port to using numpy
        points = list(points)
        def kdtree(points_left, axis=0):
            """
            kdtree(list, int) -> node

            Recursive helper function to create a kdtree from a set of points.
            Returns the root of the tree.
            """
            # reached leaf
            if not len(points_left):
                return None
            
            # each level of a kdtree must have a different axis, so we cycle
            # through axis values every level
            axis %= len(points_left[0])
        
            # Sort point list and choose median as pivot element
            points_left.sort(key=lambda point: point[axis])
            median = len(points_left) // 2
         
            node = Node(numpy.array(points_left[median]), # datum
                        kdtree(points_left[:median], axis + 1), # left
                        kdtree(points_left[median + 1:], axis + 1), # right
                        axis, # axis
                        None) # parent

            # add parents to children if necessary
            for child in (node.left,node.right):
                if child:
                    child.parent = node
            return node

        self.root = kdtree(points)
    @staticmethod
    def traverse_down(p,n):
        """
        K.traverse_down([x,y,z],Node) -> Node

        Recursive function to find where the node would be found according to
        the ordering properties for the KDTree.
        """
        if not(n.left or n.right):
            return n
        axis = n.axis
        if p[axis] <= n.datum[axis]:
            if n.left:
                return KDTree.traverse_down(p,n.left)
            else: return n
        else:
            if n.right:
                return KDTree.traverse_down(p,n.right)
            else: return n

    def nearest_neighbour(self,point):
        """
        K.nearest_neighbour([x,y,z]) -> Node

        Finds the node in the k-d tree which has the smallest euclidean
        distance stored in its datum from the point given.
        """
        point = numpy.array(point)
        euc_length_squared = lambda a,b: pow(sum(pow(a-b,2)),0.5)
        lowest_dist = float("inf")
        closest_node = None


        node = KDTree.traverse_down(point,self.root)
        visited_nodes = set() # TODO: find better way to do this
        while(node.parent):
            # add nodes to visited
            visited_nodes.add(node)

            # check if distance lower
            dist = euc_length_squared(node.datum,point)
            if dist < lowest_dist:
                closest_node = node
                lowest_dist = dist

            # look at next node
            node = node.parent
            dist = euc_length_squared(node.datum,point)
            if dist < lowest_dist:
                closest_node = node
                lowest_dist = dist

            # examine child nodes
            for child in (node.left,node.right):
                if child and (child not in visited_nodes):
                    diff = abs(child.datum[child.axis] - point[child.axis])
                    if diff < lowest_dist:
                        node = KDTree.traverse_down(point,child)
        return closest_node

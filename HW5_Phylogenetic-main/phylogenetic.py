import numpy as np
import matplotlib.pyplot as plt
import json

class DendrogramTree(object):
    def __init__(self):
        self.root = None
    
    def inorder(self):
        key_list = []
        if self.root:
            num = [0]
            self.root.inorder(num, key_list)
        return key_list, num
    
    def draw(self, thresh):
        key_list, num = self.inorder()
        max_key = max(key_list)
        
        plt.figure(figsize = (10,16))
        if self.root:
            self.root.draw(max_key)
            plt.plot([thresh, thresh], [0, num[0]])
            plt.savefig("phylogenetic_tree.png")
            plt.show()
            
    
    def add(self, key):
        if self.root:
            self.root = self.root.add(key)
        else:
            self.root = TreeNode(key)
    
    def construct_dendrogram(self, species_names, distances):
        leaf_nodes = []
        
        # add a leaf node for each species
        for species in species_names:
            leaf_nodes.append(TreeNode(species))
            
        roots = {}      # keep track of each node's root
            
        uo = UnionFindOpt(len(species_names))
        
        sorted_distances = {k: v for k, v in sorted(distances.items(), key = lambda v:v[1], reverse= True)}    # sort it by value in decreasing order
        
        species_numbers = {s:i for i,s in enumerate(species_names)}         # dict
        
        for node in leaf_nodes:
            roots[species_numbers[node.key]] = node
        
        # for each key in distances
        for pair in sorted_distances.keys():
            species = pair.split(':')    # split the species names
            
            species1 = species[0]       # get first species name
            species2 = species[1]       # get second species name
            
            # if the two species are not already in the same component
            if not uo.find(species_numbers[species1], species_numbers[species2]):
                new_node = TreeNode(sorted_distances[pair])                     # create a new node who's key is the distance between the two
                
                new_node.left = roots[uo.root(species_numbers[species1])]       # make the new node's left the root fo the left species node's root
                new_node.right = roots[uo.root(species_numbers[species2])]      # make the new node's left the root fo the left species node's root
                uo.union(species_numbers[species1], species_numbers[species2])  # union the two species nodes
                
                # make the roots of the two species node's equal to the new node
                roots[uo.root(uo.root(species_numbers[species1]))] = new_node
                roots[uo.root(uo.root(species_numbers[species2]))] = new_node
        
            
            # set the root to be the last merged node
            self.root = new_node   
    
    def get_clusters(self, thresh):
        clusters= []
        
        if self.root:
            self.root.get_clusters(thresh, clusters)
        
        return clusters

class TreeNode(object):
    def __init__(self, key):
        self.key = key
        self.left = None
        self.right = None
        self.inorder_pos = 0
    
    def add(self, key):
        ret = self
        if key < self.key:
            if self.left:
                self.left = self.left.add(key)
            else:
                self.left = TreeNode(key)
        elif key > self.key:
            if self.right:
                self.right = self.right.add(key)
            else:
                self.right = TreeNode(key)
        return ret
    
    def inorder(self, num, key_list):
        """
        Parameters
        ----------
        num: list
            List of a single element which keeps 
            track of the number I'm at
        """
        if self.left:
            self.left.inorder(num, key_list)
        self.inorder_pos = num[0]
        
        if type(self.key) is not str:
            key_list.append(self.key)
        num[0] += 1
        if self.right:
            self.right.inorder(num, key_list)
    
    def draw(self, max_key):
        y = self.inorder_pos
        
        if type(self.key) is not str:
            x = self.key
        else:
            x = max_key
            
        plt.scatter([x], [y], 50, 'k')
        plt.text(x+10, y, "{}".format(self.key))
        y_next = y-1
        
        if self.left:
            y_next = self.left.inorder_pos
            if type(self.left.key) is not str:
                x_next = self.left.key
            else:
                x_next = max_key
                
            plt.plot([x, x_next], [y, y_next])
            
            self.left.draw(max_key)
        if self.right:
            y_next = self.right.inorder_pos
            if type(self.right.key) is not str:
                x_next = self.right.key
            else:
                x_next = max_key
                
            plt.plot([x, x_next], [y, y_next])
            self.right.draw(max_key)
        
    def get_leaf_nodes(self, leaf_nodes):
        if type(self.key) is not str:
            leaf_nodes = self.left.get_leaf_nodes(leaf_nodes)
            leaf_nodes = self.right.get_leaf_nodes(leaf_nodes)
        else:
            leaf_nodes.append(self.key)
        
        return leaf_nodes
    
    def get_clusters(self, thresh, clusters):
        if type(self.key) is not str:
            if self.key < thresh:
                clusters = self.left.get_clusters(thresh, clusters)
                clusters = self.right.get_clusters(thresh, clusters)
            else:
                leaf_nodes = []
                leaf_nodes = self.get_leaf_nodes(leaf_nodes)
                clusters.append(leaf_nodes)
        else:
            clusters.append([self.key])
        
        return clusters
        

class UnionFindOpt:
    def __init__(self, N):
        self.N = N
        self.N = N
        self.parents = []
        self.weights = [1] * self.N  # create list of weights

        for i in range(N):
            self.parents.append(i)

        self._operations = 0
        self._calls = 0

    def root(self, i):
        if self.parents[i] != i:
            self._operations += 2
            self.parents[i] = self.root(self.parents[i])

        return self.parents[i]

    def union(self, i, j):
        self._calls += 1

        root_i = self.root(i)
        root_j = self.root(j)

        if root_i != root_j:
            self._operations += 1
            if self.weights[root_i] > self.weights[root_j]:
                self.parents[root_j] = root_i
                self.weights[root_i] += self.weights[root_j]
            else:
                self.parents[root_i] = root_j
                self.weights[root_j] += self.weights[root_i]

    def find(self, i, j):
        self._calls += 1
        root_i = self.root(i)
        root_j = self.root(j)

        if root_i == root_j:
            return True
        else:
            return False

    
def load_blosum(filename):
    """
    Load in a BLOSUM scoring matrix for Needleman-Wunsch

    Parameters
    ----------
    filename: string
        Path to BLOSUM file
    
    Returns
    -------
    A dictionary of {string: int}
        Key is string, value is score for that particular 
        matching/substitution/deletion
    """
    fin = open(filename)
    lines = [l for l in fin.readlines() if l[0] != "#"]
    fin.close()
    symbols = lines[0].split()
    X = [[int(x) for x in l.split()] for l in lines[1::]]
    X = np.array(X, dtype=int)
    N = X.shape[0]
    costs = {}
    for i in range(N-1):
        for j in range(i, N):
            c = X[i, j]
            if j == N-1:
                costs[symbols[i]] = c
            else:
                costs[symbols[i]+symbols[j]] = c
                costs[symbols[j]+symbols[i]] = c
    return costs

# computes the needleman wunsch distance, which can vary depending on the characters
def needleman_wunsch(str1, str2, scores):
    M = len(str1)
    N= len(str2)
    values = np.zeros((M+1, N+1))
    
    # fill in left side
    
    # BASE CASES
    # for i in the top row indices
    for i in range(1, M+1):
        # value in the next column is equal to the scores value of the next character
        # in the string plus whatever is in the cell to its left
        values[i][0] = scores[str1[i-1]] + values[i-1][0]
    
    # for i in the left column indices
    for i in range(1, N+1):
        # value in the next row is equal to the scores value of the next character
        # in the string plus whatever is in the cell above
        values[0][i] = scores[str2[i-1]]+ values[0][i-1]
        
    # INNER VALUES
    
    # for i in the row's indices
    for i in range(1, M+1):
        # for j in the column's indices
        for j in range(1, N+1):
            cost_above = values[i-1][j]         # get the cost calculated in the cell above
            cost_left = values[i][j-1]          # get the cost calculated in the cell to the left
            cost_diagonal = values[i-1][j-1]    # get the cost calculated in the cell up and to the left
            
            cost_above += scores[str1[i-1]]                 # add the cost from the scores table of the next character
            cost_left += scores[str2[j-1]]                  # add the cost from the scores table of the next character
            cost_diagonal += scores[str1[i-1]+str2[j-1]]    # add the cost from the scores table of the next character
            
            values[i][j] = max(cost_above, cost_left, cost_diagonal)    # fill the current cell with the max value of the three
            
    # print(values)
    # return the bottom right value as our needleman-wunsch value
    return values[-1, -1]

# computes the needlman wunsch values between every pair of species
def all_pairs_Needleman_Wunsch():
    species = json.load(open("organisms.json"))     # get our dictionary of {species_name: DNA Sequence}
    costs = load_blosum("blosum62.bla")             # get our dictionary of costs
    data = {}                                       # create new dict to store needleman-wunsch values in
    
    species_names = list(species.keys())            # get a list of species names to use as keys
    
    # for each species index and name
    for i, key in enumerate(species_names):
        # for each next species
        for j in range(i+1, len(species_names)):
            key2 = species_names[j]     # get the next species name
            seq1 = species[key]         # get the first species' sequence
            seq2 = species[key2]        # get the second species' sequence
            
            needleman_wunsch_distance = needleman_wunsch(seq1, seq2, costs)     # compute Needleman-Wunsch between them
            
            # store the distance between them in data dict
            data[key + ":" + key2] = needleman_wunsch_distance
    
    # dump data into json file
    json.dump({"data":data}, open("distances.json", "w"))
        
species = json.load(open("organisms.json"))     # read in the species dictionary
species_names = list(species.keys())            # get the list of species names

distances = json.load(open("distances.json"))   # read in our dictionary
distances = distances['data'] # get the only dictionary

tree = DendrogramTree()
tree.construct_dendrogram(species_names, distances)
thresh = 1260
tree.draw(thresh)
clusters = tree.get_clusters(thresh)

for element in clusters:
    print(element)

#construct_dendrogram(species_names, distances)
#dendrogram.draw()   # draw the new dendrogram
#thresh = 0.8
#dendrogram.get_clusters(thresh)
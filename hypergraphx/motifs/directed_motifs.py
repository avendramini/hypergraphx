
from hypergraphx import Hypergraph
from hypergraphx.generation.configuration_model import configuration_model
from hypergraphx.motifs.utils import (
    _directed_motifs_ho_full,
    _motifs_ho_not_full,
    _motifs_standard,
    diff_sum, norm_vector,
)

def compute_directed_motifs(hypergraph: Hypergraph, order=3, runs_config_model=10):
    """
    Compute the number of motifs of a given order in a hypergraph.

    Parameters
    ----------
    hypergraph : Hypergraph
        The hypergraph of interest
    order : int
        The order of the motifs to compute
    runs_config_model : int
        The number of runs of the configuration model

    Returns
    -------
    dict
        keys: 'observed', 'config_model', 'norm_delta'
        'observed' reports the number of occurrences of each motif in the observed hypergraph
        'config_model' reports the number of occurrences of each motif in each sample of the configuration model
        'norm_delta' reports the norm of the difference between the observed and the configuration model
        
    """

    def _motifs_order_3(edges):
        full, visited = _directed_motifs_ho_full(edges, 3)

        res = []
        for i in range(len(full)):
            res.append((full[i][0], full[i][1]))

        return res

    #def _motifs_order_4(edges):
        #full, visited = _motifs_ho_full(edges, 4)
        #not_full, visited = _motifs_ho_not_full(edges, 4, visited)
        #standard = _motifs_standard(edges, 4, visited)

        #res = []
        #for i in range(len(full)):
            #res.append((full[i][0], max([full[i][1], not_full[i][1], standard[i][1]])))

        #return res
    
    edges = hypergraph.get_edges(up_to=order)
    output = {}
    print(hypergraph)

    print("Computing observed motifs of order {}...".format(order))

    if order == 3:
        output['observed'] = _motifs_order_3(edges)
        
    elif order == 4:
        output['observed']# = _motifs_order_4(edges)
    else:
        raise ValueError("Exact computation of motifs of order > 4 is not available.")


    return output

import itertools
import math
from itertools import combinations, permutations
import time

def _motifs_ho_not_full(edges, N, visited):
    mapping, labeling = generate_motifs(N)

    T = {}
    graph = {}
    for e in edges:
        if len(e) >= N:
            continue

        T[tuple(sorted(e))] = 1

        for e_i in e:
            if e_i in graph:
                graph[e_i].append(e)
            else:
                graph[e_i] = [e]

    def count_motif(nodes):
        nodes = tuple(sorted(tuple(nodes)))
        p_nodes = power_set(nodes)

        motif = []
        for edge in p_nodes:
            if len(edge) >= 2:
                edge = tuple(sorted(list(edge)))
                if edge in T:
                    motif.append(edge)

        m = {}
        idx = 1
        for i in nodes:
            m[i] = idx
            idx += 1

        labeled_motif = []
        for e in motif:
            new_e = []
            for node in e:
                new_e.append(m[node])
            new_e = tuple(sorted(new_e))
            labeled_motif.append(new_e)
        labeled_motif = tuple(sorted(labeled_motif))

        if labeled_motif in labeling:
            labeling[labeled_motif] += 1

    for e in edges:
        if len(e) == N - 1:
            nodes = list(e)

            for n in nodes:
                for e_i in graph[n]:
                    tmp = list(nodes)
                    tmp.extend(e_i)
                    tmp = list(set(tmp))
                    if len(tmp) == N and not (tuple(sorted(tmp)) in visited):
                        visited[tuple(sorted(tmp))] = 1
                        count_motif(tmp)

    out = []

    for motif in mapping.keys():
        count = 0
        for label in mapping[motif]:
            count += labeling[label]

        out.append((motif, count))

    out = list(sorted(out))

    D = {}
    for i in range(len(out)):
        D[i] = out[i][0]

    # with open('motifs_{}.pickle'.format(N), 'wb') as handle:
    # pickle.dump(D, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return out, visited

def _directed_motifs_ho_not_full(edges, N, visited):
    mapping={}
    labeling={}
    T = {}
    graph={}
    for e in edges:
        T[tuple((tuple(sorted(e[0])),tuple(sorted(e[1]))))] = 1
        for e_i in e[0]:
            if e_i in graph:
                graph[e_i].append(e)
            else:
                graph[e_i] = [e]
        for e_i in e[1]:
            if e_i in graph:
                graph[e_i].append(e)
            else:
                graph[e_i] = [e]
            
    inverse_mapping={}
    def count_motif(nodes):
        nodes = tuple(sorted(tuple(nodes)))
        p_nodes = _all_directed_hyperedges(nodes)
        
        
        motif = []
        for edge in p_nodes:
            if edge in T:
                motif.append(edge)
        
        m = {}
        idx = 1
        for i in nodes:
            m[i] = idx
            idx += 1

        labeled_motif = []
        for e in motif:
            new_e0 = []
            for node in e[0]:
                new_e0.append(m[node])
            new_e1=[]
            for node in e[1]:
                new_e1.append(m[node])
            
            new_e = tuple((tuple(sorted(new_e0)),tuple(sorted(new_e1))))
            labeled_motif.append(new_e)
        labeled_motif = tuple(sorted(labeled_motif))
        
        if labeled_motif in labeling:
            labeling[labeled_motif] += 1
        else:
            labeling[labeled_motif]=1
        
        if labeled_motif in inverse_mapping:
            return
    

        vettore = list(range(1,N+1))
        permutazioni_vettore = permutations(vettore)
        m={}
        for permutazione in permutazioni_vettore:
            i=1
            for x in permutazione:
                m[i]=x
                i+=1
            
            new_comb=[]
            for x in labeled_motif:
                arco=[]
                for y in x:
                    parte_arco=[]
                    for j in y:
                        parte_arco.append(m[j])
                    arco.append(tuple(sorted(parte_arco)))
                
                arco=tuple(arco)
                new_comb.append(arco)
            new_comb=tuple(sorted(new_comb))
            if new_comb in mapping:
                inverse_mapping[labeled_motif]=new_comb
                mapping[new_comb].add(labeled_motif)
                return
        
        mapping[labeled_motif]=set()
        mapping[labeled_motif].add(labeled_motif)
        inverse_mapping[labeled_motif]=labeled_motif
            
        

    for e in edges:
        
        if len(e[0])+len(e[1]) == N - 1:
            nodes = list(e[0])+list(e[1])

            for n in nodes:
                for e_i in graph[n]:
                    tmp = list(nodes)
                    tmp.extend(e_i[0])
                    tmp.extend(e_i[1])
                    tmp = list(set(tmp))
                    if len(tmp) == N and not (tuple(sorted(tmp)) in visited):
                        visited[tuple(sorted(tmp))] = 1
                        count_motif(tmp)
        elif N==5 and len(e[0])+len(e[1])==N-2:
            
            nodes=list(e[0])+list(e[1])
            
            for n in nodes:
                for e_i in graph[n]:
                    if len(e_i[0])+len(e_i[1])<=N-2:
                        tmp=list(nodes)
                        tmp.extend(e_i[0])
                        tmp.extend(e_i[1])
                        tmp=list(set(tmp))
                        if len(tmp) == N and not (tuple(sorted(tmp)) in visited):
                            visited[tuple(sorted(tmp))] = 1
                            count_motif(tmp)
                        elif len(tmp)==N-1:
                            attuali=tmp.copy()
                            for second in attuali:
                                for e_second in graph[second]:
                                    if len(e_second[0])+len(e_second[1])<=N-2:
                                        tmp2=list(attuali)
                                        tmp2.extend(e_second[0])
                                        tmp2.extend(e_second[1])
                                        tmp2=list(set(tmp2))
                                        if len(tmp2)==N and not (tuple(sorted(tmp2)) in visited):
                                            visited[tuple(sorted(tmp2))] = 1
                                            count_motif(tmp2)
                                    
            
                                    
                                    
                            

    out = []

    for motif in mapping.keys():
        count = 0
        for label in mapping[motif]:
            count += labeling[label]

        out.append((motif, count))

    out = list(sorted(out))

    D = {}
    for i in range(len(out)):
        D[i] = out[i][0]

    # with open('motifs_{}.pickle'.format(N), 'wb') as handle:
    # pickle.dump(D, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return out, visited

def _motifs_standard(edges, N, visited):
    mapping, labeling = generate_motifs(N)

    graph = {}
    T = {}

    z = set()
    for e in edges:
        for n in e:
            z.add(n)

    for e in edges:
        if len(e) == 2:
            T[tuple(sorted(e))] = 1
            a, b = e
            if a in graph:
                graph[a].append(b)
            else:
                graph[a] = [b]

            if b in graph:
                graph[b].append(a)
            else:
                graph[b] = [a]

    def count_motif(nodes):
        nodes = tuple(sorted(tuple(nodes)))

        if nodes in visited:
            return

        p_nodes = power_set(nodes)

        motif = []
        for edge in p_nodes:
            edge = tuple(sorted(list(edge)))
            if edge in T:
                motif.append(edge)

        m = {}
        idx = 1
        for i in nodes:
            m[i] = idx
            idx += 1

        labeled_motif = []
        for e in motif:
            new_e = []
            for node in e:
                new_e.append(m[node])
            new_e = tuple(sorted(new_e))
            labeled_motif.append(new_e)
        labeled_motif = tuple(sorted(labeled_motif))

        if labeled_motif in labeling:
            labeling[labeled_motif] += 1

    def graph_extend(sub, ext, v, n_sub):

        if len(sub) == N:
            count_motif(sub)
            return

        while len(ext) > 0:
            w = ext.pop()
            tmp = set(ext)

            for u in graph[w]:
                if u not in sub and u not in n_sub and u > v:
                    tmp.add(u)

            new_sub = set(sub)
            new_sub.add(w)
            new_n_sub = set(n_sub).union(set(graph[w]))
            graph_extend(new_sub, tmp, v, new_n_sub)

    c = 0

    k = 0
    for v in graph.keys():
        v_ext = set()
        for u in graph[v]:
            if u > v:
                v_ext.add(u)
        k += 1
        #if k % 5 == 0:
            #print(k, len(z))

        graph_extend(set([v]), v_ext, v, set(graph[v]))
        c += 1

    out = []

    for motif in mapping.keys():
        count = 0
        for label in mapping[motif]:
            count += labeling[label]

        out.append((motif, count))

    out = list(sorted(out))

    D = {}
    for i in range(len(out)):
        D[i] = out[i][0]

    # with open('motifs_{}.pickle'.format(N), 'wb') as handle:
    # pickle.dump(D, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return out


def _motifs_ho_full(edges, N):
    mapping, labeling = generate_motifs(N)

    T = {}
    for e in edges:
        T[tuple(sorted(e))] = 1

    visited = {}

    def count_motif(nodes):
        nodes = tuple(sorted(tuple(nodes)))
        p_nodes = power_set(nodes)

        motif = []
        for edge in p_nodes:
            if len(edge) >= 2:
                edge = tuple(sorted(list(edge)))
                if edge in T:
                    motif.append(edge)

        m = {}
        idx = 1
        for i in nodes:
            m[i] = idx
            idx += 1

        labeled_motif = []
        for e in motif:
            new_e = []
            for node in e:
                new_e.append(m[node])
            new_e = tuple(sorted(new_e))
            labeled_motif.append(new_e)
        labeled_motif = tuple(sorted(labeled_motif))

        if labeled_motif in labeling:
            labeling[labeled_motif] += 1

    for e in edges:
        if len(e) == N:
            # print(e)
            visited[e] = 1
            nodes = list(e)
            count_motif(nodes)

    out = []

    for motif in mapping.keys():
        count = 0
        for label in mapping[motif]:
            count += labeling[label]

        out.append((motif, count))

    out = list(sorted(out))

    D = {}
    for i in range(len(out)):
        D[i] = out[i][0]

    # with open('motifs_{}.pickle'.format(N), 'wb') as handle:
    # pickle.dump(D, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return out, visited

def _directed_motifs_ho_full(edges, N):

    mapping={}
    labeling={}
    T = {}
    for e in edges:
        T[tuple((tuple(sorted(e[0])),tuple(sorted(e[1]))))] = 1

    visited = {}
    inverse_mapping={}
    def count_motif(nodes):
        nodes = tuple(sorted(tuple(nodes)))
        p_nodes = _all_directed_hyperedges(nodes)
        
        
        motif = []
        for edge in p_nodes:
            if edge in T:
                motif.append(edge)
        
        m = {}
        idx = 1
        for i in nodes:
            m[i] = idx
            idx += 1

        labeled_motif = []
        for e in motif:
            new_e0 = []
            for node in e[0]:
                new_e0.append(m[node])
            new_e1=[]
            for node in e[1]:
                new_e1.append(m[node])
            
            new_e = tuple((tuple(sorted(new_e0)),tuple(sorted(new_e1))))
            labeled_motif.append(new_e)
        labeled_motif = tuple(sorted(labeled_motif))
        
        if labeled_motif in labeling:
            labeling[labeled_motif] += 1
        else:
            labeling[labeled_motif]=1
        
        if labeled_motif in inverse_mapping:
            return
        
        vettore = list(range(1,N+1))
        permutazioni_vettore = permutations(vettore)
        m={}
        for permutazione in permutazioni_vettore:
            i=1
            for x in permutazione:
                m[i]=x
                i+=1
            
            new_comb=[]
            for x in labeled_motif:
                arco=[]
                for y in x:
                    parte_arco=[]
                    for j in y:
                        parte_arco.append(m[j])
                    arco.append(tuple(sorted(parte_arco)))
                
                arco=tuple(arco)
                new_comb.append(arco)
            new_comb=tuple(sorted(new_comb))
            if new_comb in mapping:
                inverse_mapping[labeled_motif]=new_comb
                mapping[new_comb].add(labeled_motif)
                return
        
        mapping[labeled_motif]=set()
        mapping[labeled_motif].add(labeled_motif)
        inverse_mapping[labeled_motif]=labeled_motif
        
        
            
                
            
        
        
        
    for e in edges:

        if len(set(e[0]+e[1])) == N and len(e[0])+len(e[1])==N and  not tuple(sorted(list(e[0])+list(e[1]))) in visited:
            nodes = list(e[0])+list(e[1])
            visited[tuple(sorted(list(e[0])+list(e[1])))]=1
            count_motif(nodes)
            
            

    out = []

    for motif in mapping.keys():
        count = 0
        for label in mapping[motif]:
            count += labeling[label]

        out.append((motif, count))

    out = list(sorted(out))

    D = {}
    for i in range(len(out)):
        D[i] = out[i][0]

    # with open('motifs_{}.pickle'.format(N), 'wb') as handle:
    # pickle.dump(D, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return out, visited

def diff_sum(observed: list, null_models: list,directed=False):
    """
    Compute the relative abundance between the observed frequencies and the null models

    Parameters
    ----------
    observed : list
        Observed frequencies
    null_models : list
        Null models

    Returns
    -------
    list
        Relative abundance between the observed frequencies and the null models

    Notes
    -----
    The relative abundance is computed as: (observed - null) / (observed + null + 4)

    """
    u_null = avg(null_models,directed)
    res = []
    if not directed:
        for i in range(len(observed)):
            res.append((observed[i][1] - u_null[i]) / (observed[i][1] + u_null[i] + 4))
    else :
        for i in range(len(observed)):
            if observed[i][0] in u_null:
                res.append((observed[i][1]-u_null[observed[i][0]])/(observed[i][1]+u_null[observed[i][0]]+4))
            else:
                res.append((observed[i][1]) / (observed[i][1] + 4))

    return res


def norm_vector(a):
    """
    Normalize a vector

    Parameters
    ----------
    a : list
        Vector to be normalized

    Returns
    -------
    list
        Normalized vector
    """
    M = 0
    for i in a:
        M += i ** 2
    M = math.sqrt(M)
    res = [i / M for i in a]
    return res


def avg(motifs, directed=False):
    
    
    if not directed:
        result = []
        for i in range(len(motifs[0])):
            s = 0
            for j in range(len(motifs)):
                #print(len(motifs[j]))
                s += motifs[j][i][1]
    
            result.append(s / len(motifs))
    else:
        m={}
        for i in range(len(motifs)):
            for j in range(len(motifs[i])):
                if motifs[i][j][0] in m:
                    m[motifs[i][j][0]]+=motifs[i][j][1]
                else:
                    m[motifs[i][j][0]]=motifs[i][j][1]
            
        result=m
    return result


def sigma(motifs):
    u = avg(motifs)

    result = []
    for i in range(len(motifs[0])):
        s = 0
        for j in range(len(motifs)):
            s += (motifs[j][i][1] - u[i]) ** 2
        s /= len(motifs)
        s = s ** 0.5

        result.append(s)
    return result


def z_score(observed, null_models):
    """
    Compute the z-score between the observed frequencies and the null models

    Parameters
    ----------
    observed : list
        Observed frequencies
    null_models : list
        Null models

    Returns
    -------
    list
        Z-score between the observed frequencies and the null models
    """
    u_null = avg(null_models)
    sigma_null = sigma(null_models)

    z_scores = []
    for i in range(len(observed)):
        z_scores.append((observed[i][1] - u_null[i]) / (sigma_null[i] + 0.01))

    return z_scores


def power_set(A):
    """
    Compute the power set of a set

    Parameters
    ----------
    A : list
        Set

    Returns
    -------
    list
        Power set of the set
    """
    subsets = []
    N = len(A)

    for mask in range(1 << N):
        subset = []

        for n in range(N):
            if ((mask >> n) & 1) == 1:
                subset.append(A[n])

        subsets.append(subset)

    return subsets

def _all_directed_hyperedges(nodi):
    """
    Compute all directed hyperedges 

    Parameters
    ----------
    A : list
        Set

    Returns
    -------
    list
        Power set of the set
    """
    iperarchi = set()

    # Genera iperarchi con nodi di partenza e di arrivo
    for lunghezza_partenza in range(1, len(nodi)):
        for nodi_partenza in combinations(nodi, lunghezza_partenza):
            for lunghezza_arrivo in range(1, len(nodi) - lunghezza_partenza + 1):
                for nodi_arrivo in combinations(set(nodi) - set(nodi_partenza), lunghezza_arrivo):
                    iperarco = (tuple(sorted(nodi_partenza)), tuple(sorted(nodi_arrivo)))
                    iperarchi.add(iperarco)

    return iperarchi
    

def _is_connected(edges, N):
    nodes = set()
    for e in edges:
        for n in e:
            nodes.add(n)

    if len(nodes) != N:
        return False

    visited = {}
    for i in nodes:
        visited[i] = False
    graph = {}
    for i in nodes:
        graph[i] = []

    for edge in edges:
        for i in range(len(edge)):
            for j in range(len(edge)):
                if edge[i] != edge[j]:
                    graph[edge[i]].append(edge[j])
                    graph[edge[j]].append(edge[i])

    q = []
    nodes = list(nodes)
    q.append(nodes[0])
    while len(q) != 0:
        v = q.pop(len(q) - 1)
        if not visited[v]:
            visited[v] = True
            for i in graph[v]:
                q.append(i)
    conn = True
    for i in nodes:
        if not visited[i]:
            conn = False
            break
    return conn


def relabel(edges: list, relabeling: dict):
    """
    Relabel the vertices of a hypergraph according to a given relabeling

    Parameters
    ----------
    edges : list
        Edges of the hypergraph
    relabeling : dict
        Relabeling

    Returns
    -------
    list
        Edges of the hypergraph with the vertices relabeled

    Notes
    -----
    The relabeling is a dictionary that maps the old labels to the new labels
    """
    res = []
    for edge in edges:
        new_edge = []
        for v in edge:
            new_edge.append(relabeling[v - 1])
        res.append(tuple(sorted(new_edge)))
    return sorted(res)


def generate_motifs(N):
    """
    Generates all possible patterns of non-isomorphic subhypergraphs of size N

    Parameters
    ----------
    N : int
        Size of the subhypergraphs

    Returns
    -------
    list
        List of all possible patterns of non-isomorphic subhypergraphs of size N
    """
    n = N
    
    assert n >= 2

    h = [i for i in range(1, n + 1)]
    A = []

    for r in range(n, 1, -1):
        A.extend(list(itertools.combinations(h, r)))

    B = power_set(A)

    C = []
    for i in range(len(B)):
        if _is_connected(B[i], N):
            C.append(B[i])

    isom_classes = {}

    for i in C:
        edges = sorted(i)
        relabeling_list = list(itertools.permutations([j for j in range(1, n + 1)]))
        found = False
        for relabeling in relabeling_list:
            relabeling_i = relabel(edges, relabeling)
            # print(relabeling_i)
            if tuple(relabeling_i) in isom_classes:
                found = True
                break
        if not found:
            isom_classes[tuple(edges)] = 1

    mapping = {}
    labeling = {}

    for k in isom_classes.keys():
        mapping[k] = set()
        relabeling_list = list(itertools.permutations([j for j in range(1, n + 1)]))
        for relabeling in relabeling_list:
            relabeling_i = relabel(k, relabeling)
            labeling[tuple(sorted(relabeling_i))] = 0
            mapping[k].add(tuple(sorted(relabeling_i)))

    return mapping, labeling


def count_directed_iso_classes(N):
    edges=_all_directed_hyperedges(list(range(1,N+1)))
    iso_classes={}
    recursion(0,[],N,list(edges),iso_classes)
    return len(iso_classes)


def recursion(position, motif,N,edges, iso_classes):
    if position==len(edges):
        #check if connected
        
        parent=list(range(N))

        def find(x):
            if parent[x] == x:
                return x
            parent[x] = find(parent[x])  # Path compression
            return parent[x]

        def union(x, y):
            rootX = find(x)
            rootY = find(y)
            if rootX != rootY:
                parent[rootY] = rootX

        for edge in motif:
            for node in edge[0]:
                union(node-1,edge[0][0]-1)
            for node in edge[1]:
                union(node-1,edge[0][0]-1)
        
        unique_roots = set(find(x-1) for x in range(N))
        if len(unique_roots) !=1:
            return
        


        #check isomorphism class
        motif=tuple(sorted(motif))
        vettore = list(range(1,N+1))
        permutazioni_vettore = permutations(vettore)
        m={}
        trovato=False
        for permutazione in permutazioni_vettore:
            i=1
            for x in permutazione:
                m[i]=x
                i+=1
            
            new_comb=[]
            for x in motif:
                arco=[]
                for y in x:
                    parte_arco=[]
                    for j in y:
                        parte_arco.append(m[j])
                    arco.append(tuple(sorted(parte_arco)))
                
                arco=tuple(arco)
                new_comb.append(arco)
            new_comb=tuple(sorted(new_comb))
            if new_comb in iso_classes:
                trovato=True
                
        if not trovato:
            iso_classes[motif]=1
        return
    
    recursion(position+1,motif.copy(),N,edges,iso_classes)
    motif.append(edges[position])
    recursion(position+1,motif.copy(),N,edges,iso_classes)
    motif.pop()


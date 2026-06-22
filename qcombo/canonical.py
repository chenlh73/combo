# NOTE：
# use set to store up or lo indices
# use list to store all indices.
# because list can keep the order of ele while set may change the original order. 

# PROBLEMS:
# [√] The order of set elements may change irregularly and not preserve the original order? Why?
# The 'in' operator between FinishSet collections seems to have issues?


#use tuple to store indices and use list to process indices.DO NOT USEset!!

from sympy import IndexedBase,symbols
from .tools import ProgressBar
from sympy.core.add import Add



# Worker function for parallel canonicalize (must be at module level for pickling)
def _worker_canonical_order(term):
    """Worker function for ProcessPoolExecutor in canonicalize."""
    X = IndexedBase('X')
    if type(term) != type(X[{},{}]):
        coef = term.args[0]
        if type(coef) != type(X[{},{}]):
            # has coefficient
            return coef * canonicalOrder(term / coef)
        else:
            return canonicalOrder(term)
    else:
        return term


X=IndexedBase('X')
delta=IndexedBase(chr(948))
n=IndexedBase('n')


def _selectDeltaTerms(_inputTerms):
    """
    Separate delta tensors from other terms in a Mul expression.
    
    Returns:
        (_deltaTerms, _otherTerms): tuple of delta product and other product
    """
    _deltaTerms=1
    _otherTerms=1
    for i in _inputTerms.args:
        if i.base==delta:
            _deltaTerms*=i
        else:
            _otherTerms*=i
    return _deltaTerms,_otherTerms


def _formDeltaRule(deltaTerms)->list:
    """
    Generate index replacement rules from delta terms.
    
    For each delta^a_b, creates a rule to replace the larger index with the smaller.
    
    Returns:
        List of [old_index, new_index] replacement rules
    """
    if deltaTerms==1:
        return []
    res=[]
    if deltaTerms.args[0]==delta:# only one delta item 
        eleDelta=tuple([deltaTerms])
    else:
        eleDelta=deltaTerms.args#separate delta terms  of 'deltaTerms'

    for i in eleDelta:#TODO: optimize this structure
        up,lo=i.args[1],i.args[2]
        if (str(up)>str(lo)):# must trans to str to compare.
            # beenReplace=symbols(str(up))
            # afterReplace=symbols(str(lo)) #big index to small index. Is it necessary to trans to a small one?
            beenReplace = up
            afterReplace = lo
        else:
            # beenReplace=symbols(str(lo))
            # afterReplace=symbols(str(up))
            beenReplace = lo
            afterReplace = up
        res.append([beenReplace,afterReplace])
    return res

def _intersection(set1,set2):
    '''
    if set1( or list1) have common elements with set2 return true.
    '''
    for i in set1:
        for j in set2:
            if i==j: return True
    return False

def changeIndices(exp,allrules):
    '''
    Replace the indices  by appling the rule.

    ## Case:
    rule={{a->b},{c->d}}
    rule={{codition,a->b},{codition,c->d}}
    rule=[[{a},{b}],[{c,f},{d}]]
    '''
    if allrules==[]:#no need to apply rule
        return exp

    for rule in allrules:#replace rule[0] to rule[1] for each allrules' element.
        if type(exp)==type(X[(),()]):
            expEle=(exp,)
        else:
            expEle=exp.args#separate the ele of exp.
        tmp=1
        for i in expEle:
            newUp,newLo=[],[]#stor up and low indices.
            eleArgs=i.args#separate the indices of ele.
            oldUp,oldLo=eleArgs[1],eleArgs[2]
            if _intersection(rule[0],oldUp): #replace the up symbols of one element.
                for j in oldUp:
                    if j in rule[0]:
                        newUp.append(rule[1][0])
                    else: newUp.append(j)
            else:newUp=oldUp
                # for k in rule[0]:
                #  for j in oldUp:
                #     if j==k:
                #          newUp.append(rule[1])
                #     else: newUp.append(j)
            if _intersection(rule[0],oldLo):#replace the dow symbols of one element.
                for j in oldLo:
                    if j in rule[0]:
                        newLo.append(rule[1][0])
                    else: 
                        newLo.append(j)
            else: newLo=oldLo
            tmp*=eleArgs[0][tuple(newUp),tuple(newLo)]
        exp=tmp

    return exp

def processDeltaTerms(inputTerms):
    '''
    Process the delta terms of the expression.

    # Examples:
    >>> processDeltaTerms(G[{i,j},{a,b}]*H[{c,d},{k,l}]*delta[{a},{c}]*delta[{b},{d}])
    G[{i, j}, {a, b}] H[{a, b}, {k, l}]

    # Note:
    - When the sets have common eles, it just show one of them, {1,1}={1}.
    - Will the inputTerms' indices have common ele?or appear after i process the delta  terms?
    - notice that delta[left,right]'s left and right come from idxup and idxLo,so it is impossible to 
    appear the common ele in a up or lo indices sets. It means that H[{c,c},{k,l}] can not appear.

    '''
    deltaTerms,otherTerms=_selectDeltaTerms(inputTerms)
    deltaRule=_formDeltaRule(deltaTerms)
    deltaRemoved=changeIndices(otherTerms,deltaRule)
    return deltaRemoved


def processDelta():
    '''
    Prcess delta terms of the exp when exp is the sum of multiterms like: A{{},{}}delta[{},{}]+B[{},{}]delta[{},{}]
    '''
    pass

def _getIndices(exp):
    res=[]
    if type(exp)==type(X[(),()]):
        res=[[exp.args[1],exp.args[2]]]
    else:
        for i in exp.args:
            up=i.args[1]
            lo=i.args[2]
            res+=[[up,lo]]
    return res

def _flattenIdx(indices):
    """
    Flatten nested index structure into a single list.
    
    Example: [[(a,b), (c,d)], [(e,), (f,)]] -> [a, b, c, d, e, f]
    """
    res=[]
    for i in indices:
        for j in i:
            res+=list(j)
    return res

def _internalAndRepeatedIdx(flatIdx,externalIdx=[]):
    """
    Find indices that appear exactly twice (dummy/internal indices).
    
    Parameters:
        flatIdx: flattened list of all indices
        externalIdx: list of external (free) indices to exclude
    
    Returns:
        List of repeated (dummy) indices
    """
    res=[]
    idxdict={}#initial a dict
    for i in flatIdx:
        idxdict[i]=0
    for i in flatIdx:
        if i not in externalIdx:
            idxdict[i]+=1
            if(idxdict[i]==2):
                res.append(i)
    return res

def _getindexSet(indices,internalIdx):
    """
    Extract index subsets that contain internal (dummy) indices.
    
    Returns:
        List of index tuples containing internal indices
    """
    indexSets=[]
    for i in indices:
        for j in i:
            if(_union(j,internalIdx)):
                indexSets.append(_union(j,internalIdx)) 
    return indexSets


def _union(set1,set2):
    """
    Find intersection of two sequences (common elements).
    
    Returns:
        Tuple of common elements
    """
    res=[]
    for i in set1:
        for j in set2:
            if i ==j:
                res.append(i)
    return tuple(res)#it seems that the set cannot keep the original order.#yes so i changed it into tuple

def _diff(set1,set2):
    '''return the diffent ele of set2 with set1'''
    res=[]
    for i in set2:
        tmp=0
        for j in set1:
            if j==i:
                tmp=1
        if tmp==0:
            res.append(i)
    return tuple(res)    

def _findEquivalentIndices(indexSets):
    """
    Find groups of equivalent indices that can be interchanged.
    
    Removes common elements between index sets to identify independent groups.
    
    Returns:
        List of index groups that can be independently renamed
    """
    equi=[]
    selectedSet=set()#stored selected ele in one-dimension set.
    for beingSelectedSet in indexSets:
        diffset=_diff(selectedSet,beingSelectedSet)
        if diffset!=():
            equi.append(diffset)
            selectedSet.update(set(diffset))# set will remove the common elements when you update it.
    return equi

def _replaceRules(equivalents):
    """
    Generate replacement rules mapping equivalent indices to dummy variables.
    
    Assigns alphabetical dummy names (aa, bb, cc, ...) to equivalent index groups.
    
    Returns:
        List of replacement rules [[old_index, new_index], ...]
    """
    dummyPos=97 # chr(97) = 'a'
    rules=[]
    for set in equivalents:
        for j in set:
            rules.append([(j,),(symbols(chr(dummyPos)*2),)])
            dummyPos+=1
    return rules

def processEmptySetIn_n(exp):
    '''
    Replace EmptySet indices in 'n' to {}.
    '''
    res=1
    if type(exp)==type(X[(),()]):
        if exp.base==n:
            res*=n[(),exp.args[2]]
        else:res=exp
    else:
        for i in exp.args:
            if i.base==n:
                res*=n[(),i.args[2]]
            else:res*=i     
    return res

def canonicalOrder(term):
    """
    Canonicalize a single term by standardizing dummy index names.
    
    Steps:
        1. Process delta terms (contract indices)
        2. Find internal (dummy) indices
        3. Generate replacement rules
        4. Apply rules to rename dummy indices
        5. Fix empty index sets in occupation numbers
    
    Returns:
        Term with standardized dummy indices
    """
    res=processDeltaTerms(term)# processdeltaTerms terms
    indices=_getIndices(res)# extract out  all indices from res.
    flatIdx= _flattenIdx(indices)# all of indices in a line.
    internalIdx=_internalAndRepeatedIdx(flatIdx,externalIdx=[])#find reaped elements in flatIdx.
    indexSets=_getindexSet(indices,internalIdx)#find out the reaped indices and keep original order.
    equivalents=_findEquivalentIndices(indexSets)#remove common elements from the reaped indices which keep the original order 
    rules=_replaceRules(equivalents)#form the replace fules.
    res=changeIndices(res,rules)#apply rules
    res=processEmptySetIn_n(res)
    return res





# ========== Original canonicalize (before parallel optimization) ==========
def canonicalize_old(exp):
    res=0
    progress = ProgressBar(len(exp.args), "canonicalizing")
    for i in  exp.args:
        if type(i)!=type(X[{},{}]):
            coef=i.args[0]
            if type(coef)!=type(X[{},{}]):#if it has coefficient,
                res+=coef*canonicalOrder(i/coef)# get rid of the coefficient before use canonicalOrder function.
            else:
                res+=canonicalOrder(i)
        else:
            res+=i
        progress.update()
    return res
# ========================================================================


def canonicalize(exp, parallel=False, show_process=True):
    """
    Canonicalize expression by renaming dummy indices.
    
    Args:
        exp: SymPy expression (Add)
        parallel: Enable parallel computing using ProcessPoolExecutor
        show_process: Show progress bar
    """
    terms = list(exp.args)
    total = len(terms)
    
    if parallel and total > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        # Parallel computing
        if show_process:
            print(f"Parallel canonicalizing {total} terms...")
        
        results = []
        with ProcessPoolExecutor() as executor:
            # Submit all tasks
            future_to_idx = {executor.submit(_worker_canonical_order, term): i 
                             for i, term in enumerate(terms)}
            
            # Collect results with progress bar
            if show_process:
                progress = ProgressBar(total, "parallel canonicalizing")
            
            for future in as_completed(future_to_idx):
                results.append(future.result())
                
                if show_process:
                    progress.update()
        
        # Sum up all results
        res = Add(*results)
        return res
    else:
        # Sequential computing (original implementation)
        res = 0
        if show_process:
            progress = ProgressBar(total, "canonicalizing")
        
        for i in terms:
            if type(i) != type(X[{},{}]):
                coef = i.args[0]
                if type(coef) != type(X[{},{}]):
                    res += coef * canonicalOrder(i / coef)
                else:
                    res += canonicalOrder(i)
            else:
                res += i
            if show_process:
                progress.update()
        
        return res



if __name__=='__main__':
    from sympy import IndexedBase
    from sympy import symbols

    Alpha=IndexedBase(chr(913))
    xi=IndexedBase(chr(958))
    lamda=IndexedBase(chr(955))
    delta=IndexedBase(chr(948))
    n=IndexedBase('n')

    X=IndexedBase('X')
    A=IndexedBase('A')
    B=IndexedBase('B')
    C=IndexedBase('C')
    G=IndexedBase('G')
    H=IndexedBase('H')
    a=IndexedBase('a')
    b=IndexedBase('b')
    p=IndexedBase('p')
    q=IndexedBase('q')
    r=IndexedBase('r')
    s=IndexedBase('s')
    

    case1=G[(a,), (b,)]*H[(p, q), (r, s)]*delta[(a,), (s,)]*lamda[(p, q), (b, r)]
    case1=A[(p,), (s,)]*G[(a,), (b,)]*H[(p, q), (r, s)]*n[(), (q,)]*delta[(a,), (r,)]*delta[(q,), (b,)]
    case1=G[(a,), (b,)]*H[(p, q), (r, s)]*A[(p, q), (b, s)]*delta[(a,), (r,)]
    # case1=H[(p, q), (r, s)]*delta[(a,), (s,)]
    #G[{aa}, {bb}]*H[{cc, dd}, {aa, ee}]*λ[{cc, dd}, {bb, ee}]
    caseRes=canonicalOrder(case1)
    print(case1)
    print(caseRes)
    # inputcase=A[{i, j}, {a, b}]*B[{a, b}, {k, l}] + 5*A[{i, j}, {c, d}]*B[{c, d}, {k, l}] + lamda[{c, d}, {e, f}]*A[{i, j}, {c, d}]*B[{e, f}, {k, l}]+A[{1},{2}]
    


if __name__=='canonical':
    from sympy import IndexedBase
    from sympy import symbols
    X=IndexedBase('X')
    delta=IndexedBase(chr(948))
    n=IndexedBase('n')
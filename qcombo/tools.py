from sympy import IndexedBase #for storing indices
from sympy import simplify,expand
from sympy import symbols,preorder_traversal,factorial
from sympy import S
from sympy.tensor.indexed import Indexed
from sympy.core.mul import Mul
from sympy.core.add import Add
import time
from itertools import permutations
import qcombo




#show progress  update by Chen.L.H 2025/12/22 
class ProgressBar:
    """
    Simple text-based progress bar with time estimation.
    
    Example:
        progress = ProgressBar(100, "Processing")
        for i in range(100):
            # do work
            progress.update()
    """
    def __init__(self, total, description="Doing...", bar_length=50):
        self.total = total
        self.description = description
        self.bar_length = bar_length
        self.current = 0
        self.start_time = time.time()
        
    def update(self):
        self.current += 1
        progress = self.current / self.total
        progress_info = f"[{self.current}/{self.total}]{progress:.1%}"
        
        # Progress bar
        filled_length = int(self.bar_length * progress)
        bar = '█' * filled_length + '░' * (self.bar_length - filled_length)
        
        # Time estimation
        elapsed = time.time() - self.start_time
        if self.current > 1:
            estimated_total = elapsed / progress
            remaining = estimated_total - elapsed
            time_info = f" | Remaining: {remaining:.1f}s"
        else:
            time_info = ""
        
        # Use \r to return to line start, and add flush=True to force immediate output
        print(f'\r{self.description}: [{bar}] {progress_info} {time_info}', end='', flush=True)
        

        if self.current == self.total:
            print(f'\n{self.description} completed! Total Time:{elapsed:.1f}s')



X=IndexedBase('X')
Y=IndexedBase('Y')#  X,Y just used to present the type of different input expression.like some multi-expr:X[{},{}]*Y[{},{}] or add-expr:X[{},{}]+Y[{},{}] 
A=IndexedBase('A')# A used to identify the exist of it in expression.
xi=IndexedBase(chr(958))
lamda=IndexedBase(chr(955))
delta=IndexedBase(chr(948))
n=IndexedBase('n')


class SimplifyRule:
    """
    Simplification rules for Wick theorem expressions.
    
    Provides methods to convert xi (contraction) tensors to lambda-delta form,
    and to diagonalize single-particle density matrices.
    """

    # Define the symbols used in the rules
    # xi=IndexedBase(chr(958))
    # lamda=IndexedBase(chr(955))
    # delta=IndexedBase(chr(948))
    # n=IndexedBase('n')

    def xiRule_oold(exp, show_process=True):
        """
        Replace xi tensors with (lambda - delta). [Legacy version]
        
        xi^a_b = lambda^a_b - delta^a_b
        """
        expArgs=exp.args#split expression by addition into tuple elements
        res=0
        if show_process:
            progress = ProgressBar(len(expArgs), "xiRule applying")
        for i in expArgs:
            elemArgs=i.args #split each element by multiplication into sub-elements
            subexp=1
            for j in elemArgs:
                if (j!=-1 and j.base==xi):
                    superscript=j.args[1]
                    subscript=j.args[2]
                    j=lamda[superscript,subscript]-delta[superscript,subscript]
                subexp*=j
            res+=subexp.expand()
            if show_process:
                progress.update()
        
        return simplify(res)


    def natRule_oold(exp, show_process=True):
        """
        Diagonalize single-particle density matrix. [Legacy version]
        
        lambda^a_b (single-body) = n_a * delta^a_b
        """
        expArgs=exp.args#split expression by addition into tuple elements
        res=0 
        if show_process:
            progress = ProgressBar(len(expArgs), "natRule applying")
        for i in expArgs:
            elemArgs=i.args#split each element by multiplication into sub-elements
            subexp=1
            for j in elemArgs:
                if (j!=-1 and j.base==lamda and len(j.args[1])==1):
                    superscript=j.args[1]
                    subscript=j.args[2]
                    j=n[(),superscript]*delta[superscript,subscript]
                subexp*=j
            res+=subexp
            if show_process:
                progress.update()
                
        return simplify(res)

    
    # New xiRule and natRule using xreplace
    # Compared to the original method that iterates each term and searches for lambda/xi,
    # xreplace with preorder_traversal is faster.

    def xiRule_old(expr, show_process=True):
        """
        Replace xi tensors with (lambda - delta) using preorder traversal.
        
        xi^a_b = lambda^a_b - delta^a_b
        
        Parameters:
            expr: SymPy expression
            show_process: Show progress bar
        
        Returns:
            Simplified expression with xi tensors replaced
        """
        replacements = {}
        if show_process:
            progress = ProgressBar(len(list(preorder_traversal(expr))), "xiRule applying")
        for term in preorder_traversal(expr):
            if isinstance(term, Indexed) and term.base == xi:
                up_idx, lo_idx = term.indices
                replacements[term] = (lamda[up_idx, lo_idx] - delta[up_idx, lo_idx])
            if show_process:
                progress.update()
        res = expr.xreplace(replacements)
        simplify_res = expand(res)
        
        return simplify_res

    def natRule_old(expr, show_process=True):
        """
        Diagonalize single-particle density matrix using preorder traversal.
        
        lambda^a_b (single-body) = n_a * delta^a_b
        
        Parameters:
            expr: SymPy expression
            show_process: Show progress bar
        
        Returns:
            Expression with single-body lambda replaced by occupation numbers
        """
        replacements = {}
        if show_process:
            progress = ProgressBar(len(list(preorder_traversal(expr))), "natRule applying")
        for term in preorder_traversal(expr):
            if isinstance(term, Indexed) and term.base == lamda and len(term.args[1])==1:
                up_idx, lo_idx = term.indices
                replacements[term] = n[(),up_idx]*delta[up_idx,lo_idx]
            if show_process:
                progress.update()
        res = expr.xreplace(replacements)
        simplify_res = expand(res)
        
        return simplify_res

    # New filtering approach: no need to iterate and search for single-body lambda or xi
    # Directly iterate over all possible index combinations and replace
    # xreplace skips silently when no matching single-body lambda or xi is found
    #
    def xiRule(expr, all_indices, show_process=True):
        """
        Replace xi tensors with (lambda - delta) using index enumeration.
        
        Faster than iterative traversal method.
        
        Parameters:
            expr: SymPy expression
            all_indices: list of all possible indices to enumerate
            show_process: Show progress messages
        
        Returns:
            Simplified expression with xi tensors replaced
        """
        replacements = {}
        if show_process:
            print("xiRule applying")
        for up_idx in all_indices:
            for lo_idx in all_indices:
                replacements[xi[(up_idx,), (lo_idx,)]] = (lamda[(up_idx,), (lo_idx,)] - delta[(up_idx,), (lo_idx,)])
        res = expr.xreplace(replacements)
        simplify_res = expand(res)
        if show_process:
            print("xiRule completed!")
        return simplify_res

    def natRule(expr, all_indices, show_process=True):
        """
        Diagonalize single-particle density matrix using index enumeration.
        
        Faster than iterative traversal method.
        
        Parameters:
            expr: SymPy expression
            all_indices: list of all possible indices to enumerate
            show_process: Show progress messages
        
        Returns:
            Expression with single-body lambda replaced by occupation numbers
        """
        replacements = {}
        if show_process:
            print("natRule applying")
        for up_idx in all_indices:
            for lo_idx in all_indices:
                replacements[lamda[(up_idx,), (lo_idx,)]] = n[(),(up_idx,)]*delta[(up_idx,),(lo_idx,)]
        res = expr.xreplace(replacements)
        simplify_res = expand(res)
        if show_process:
            print("natRule completed!")
        
        return simplify_res


# Worker function for parallel filterbody
def _worker_filterbody(args):
    """Worker function for ProcessPoolExecutor in filterbody_add."""
    term, body_type = args
    if isinstance(term, Mul):
        return Filter.filterbody_mul(term, body_type)
    elif isinstance(term, Indexed):
        return Filter.filterbody_tensor(term, body_type)
    else:
        if body_type == 0:
            return term
        else:
            return 0


class Filter:

    # Define the symbols used in the filter
    X=IndexedBase('X')
    Y=IndexedBase('Y')#  X,Y just used to present the type of different input expression.like some multi-expr:X[{},{}]*Y[{},{}] or add-expr:X[{},{}]+Y[{},{}] 
    A=IndexedBase('A')# A used to identify the exist of it in expression.

    def filterbody_old(terms,bodyType, show_process=True):
        if (type(terms)!=type(X[{},{}]+Y[{},{}])):#input terms like A[{},{}]*B[{},{}] or A[{},{}]
            if(type(terms)==type(X[{},{}]*Y[{},{}])):#like A[{},{}]*B[{},{}] or -A[{},{}]*B[{},{}]
                if (type(terms.args[0])==type(A[{},{}])):#if input terms like A[{},{}]*B[{},{}]
                    firstTerm=terms.args[0]
                    if(firstTerm.base==A and len(firstTerm.args[1])==bodyType ):
                        return terms
                    else: return 0
                elif(type(terms.args[1])==type(A[{},{}])):#if input terms like -A[{},{}]*B[{},{}]
                    secondTerm=terms.args[1]
                    if(secondTerm.base==A and len(secondTerm.args[1])==bodyType ):
                        return terms
                    else: return 0    
                elif(bodyType==0):# 0-body term
                    return terms
                else: return 0 
            else:# like  A[{},{}]
                # if(terms.base==A and len(terms.args[1]==bodyType)):
                if(terms.base==A and len(terms.args[1])==bodyType):
                    return terms
                elif(bodyType==0):
                    return terms 
                else: 
                    return 0

        else:#terms like A[{},{}]*B[{},{}] + C[{},{}]*D[{},{}]
            res=0
            if show_process:
                progress = ProgressBar(len(terms.args), "filtering")
            for i in terms.args:# i:(A[{},{}]*B[{},{}] , C[{},{}],-D[{},{}]*E[{},{}])
                if(type(i)==type(X[{},{}]*Y[{},{}])):#i like A[{},{}]*B[{},{}] or -A[{},{}]*B[{},{}]
                    if (type(i.args[0])==type(A[{},{}])):#if input terms like A[{},{}]*B[{},{}]
                        firstTerm=i.args[0]
                        if(firstTerm.base==A and len(firstTerm.args[1])==bodyType ):
                            res+=i
                        elif(firstTerm.base!=A and bodyType==0):
                            res+=i
                    elif(type(i.args[1])==type(A[{},{}])):#if input terms like -A[{},{}]*B[{},{}]
                        secondTerm=i.args[1]
                        if(secondTerm.base==A and len(secondTerm.args[1])==bodyType ):
                            res+=i
                        elif(secondTerm.base!=A  and bodyType==0):
                            res+=i
                    elif(bodyType==0):
                        res+=i
                else:# like  A[{},{}]
                    if(i.base==A and len(i.args[1])==bodyType):
                        res+=i
                    elif(i.base!=A and bodyType==0):
                        res+=i
                if show_process:
                    progress.update()
            return res

    def filterbody_tensor(expr,body_type):

        if not isinstance(expr, Indexed):
            raise ValueError("Input must be an Indexed object")

        # find_A = False
        tensor_body = 0

        if expr.base == A:
            # find_A = True
            tensor_body = len(expr.indices[0])
        
        if tensor_body == body_type:
            return expr
        else:
            return 0

    def filterbody_mul(expr,body_type):
        if not isinstance(expr, Mul):
            raise ValueError("Input must be a Mul object")
        
        # find_A = False
        expr_body = 0

        for term in expr.args:
            if isinstance(term, Indexed) and term.base == A:
                # find_A = True
                expr_body = len(term.indices[0])
                break
        
        if expr_body == body_type:
            return expr
        else:
            return 0

    def filterbody_add(expr, body_type, show_process=True, parallel=False):
        """
        Filter body_type from Add expression.
        
        Args:
            expr: SymPy Add expression
            body_type: Target body type to filter
            show_process: Show progress bar
            parallel: Enable parallel computing using ProcessPoolExecutor
        """
        if not isinstance(expr, Add):
            raise ValueError("Input must be an Add object")
        
        terms = list(expr.args)
        total = len(terms)
        
        if parallel and total > 1:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            
            if show_process:
                print(f"Parallel filtering {total} terms...")
            
            # Prepare task arguments
            task_args = [(term, body_type) for term in terms]
            
            results = []
            with ProcessPoolExecutor() as executor:
                # Submit all tasks
                future_to_idx = {executor.submit(_worker_filterbody, args): i 
                                 for i, args in enumerate(task_args)}
                
                # Collect results with progress bar
                if show_process:
                    progress = ProgressBar(total, "parallel filtering")
                
                for future in as_completed(future_to_idx):
                    results.append(future.result())
                    if show_process:
                        progress.update()
            
            # Filter out None or 0 results and create Add
            filtered_results = [r for r in results if r != 0]
            return Add(*filtered_results)
        else:
            # Sequential computing (original implementation)
            res = []
            if show_process:
                progress = ProgressBar(total, "filtering")
            
            for term in terms:
                if isinstance(term, Mul):
                    res.append(Filter.filterbody_mul(term, body_type))
                elif isinstance(term, Indexed):
                    res.append(Filter.filterbody_tensor(term, body_type))
                else:
                    if body_type == 0:
                        res.append(term)

                if show_process:
                    progress.update()

            return Add(*res)

    def filterbody(expr, body_type, show_process=True, parallel=False):
        """
        Filter input body_type of Normal ordered symbol "A" from expr.
        
        Args:
            expr: SymPy expression
            body_type: Target body type to filter
            show_process: Show progress bar
            parallel: Enable parallel computing using ProcessPoolExecutor
        """

        if isinstance(expr, Add):
            return Filter.filterbody_add(expr, body_type, show_process, parallel)
        elif isinstance(expr, Mul):
            return Filter.filterbody_mul(expr, body_type)
        elif isinstance(expr, Indexed):
            return Filter.filterbody_tensor(expr, body_type)
        else:
            if body_type == 0:
                return expr
            else:
                return 0



##############################################################################

def _sepatate(tmpCom,otherTerms):
    '''
    Select terms from otherTerms which has tmpCom term  to removeComTerms.
    Store term to other terms to otherTerms.
    '''
    otherComTerms=0
    if type(otherTerms)==type(X[{},{}]+Y[{},{}]):
        otherTerms=otherTerms.args
    else: otherTerms=set([otherTerms])

    for i  in otherTerms:#(-C*D,E*F)
        for j  in i.args:#(-1,C,D)
            if (j==tmpCom):
                otherComTerms+=i/tmpCom
                break
    return otherComTerms

def uniteSimilarTerms(exp):
    """
    Combine terms with identical tensor structure.
    
    Groups Mul terms that share common tensor factors and combines
    their remaining coefficients.
    
    Parameters:
        exp: SymPy Add expression
    
    Returns:
        Expression with like terms factored and combined
    
    Example:
        A*B*X + A*B*Y + C*D*Z >> A*B*(X+Y) + C*D*Z
    """
    if(type(exp)!=type(X[{},{}]+Y[{},{}])):# return exp when it don't have add terms
        return exp
    baseTerms=exp.args[0]# select A*B from A*B+C*D
    otherTerms=exp-baseTerms
    for tmpCom in baseTerms.args:#turn A*B to (A,B),term=A,B
        com=tmpCom
        comTermsDivideCom=baseTerms/com
        otherComTermsDivideCom=_sepatate(tmpCom,otherTerms)

        if(otherComTermsDivideCom!=0):
            comTermsDivideCom+=otherComTermsDivideCom
            otherTerms=exp-(comTermsDivideCom*com).expand()
            break
    res=com*(uniteSimilarTerms(comTermsDivideCom))+uniteSimilarTerms(otherTerms)#recurse it. 
    return res

################################################################################

def _tupleMultTosimp(tupleUp,tupleLo):
    '''
    Remove the repetitive index of element in Up and Lo.
    '''
    up=[]
    lo=[]
    for i in tupleUp:
        up.append(symbols(str(i)[0]))
        indicesSet.add(str(i)[0])
    for i in tupleLo:
        lo.append(symbols(str(i)[0]))
        indicesSet.add(str(i)[0])
    return tuple(up),tuple(lo)


# This function may also need parallel processing
def indicesMultToSimp_old(exp):
    ''''
    The indices of canonicalize result are repetitive. This function is designed to simplify the indices.  
    '''
    global indicesSet
    indicesSet=set()
    if (type(exp)!=type(X[{},{}]+Y[{},{}])):#input exp like A[{},{}]*B[{},{}] or A[{},{}]
        res=1
        if(type(exp)==type(X[{},{}]*Y[{},{}])):#like A[{},{}]*B[{},{}] or -A[{},{}]*B[{},{}]
            for i in exp.args:
                if(type(i)==type(A[{},{}])):
                    res*=i.base[_tupleMultTosimp(i.args[1],i.args[2])]
                else: res*=i
        else:# like  A[{},{}] or -A[{},{}] #Impossibly occur because commute Terms would never have such terms
            i=exp
            res=i.base[_tupleMultTosimp(i.args[1],i.args[2])]
    else:#exp like A[{},{}]*B[{},{}] + C[{},{}]*D[{},{}]-E[]F[]
        res=0
        for i in exp.args:#i:(A[{},{}]*B[{},{}] ,   C[{},{}],   -D[{},{}]*E[{},{}])
            if(type(i)==type(X[{},{}]*Y[{},{}])):#i like A[{},{}]*B[{},{}] or -A[{},{}]*B[{},{}]
                tmp=1
                for j in i.args:#j (-1,A,B)
                    if(type(j)==type(A[{},{}])):
                        tmp*=j.base[_tupleMultTosimp(j.args[1],j.args[2])]
                    else:tmp*=j
                res+=tmp
            else:# likeA[{},{}]
                res+=i.base[_tupleMultTosimp(i.args[1],i.args[2])]
    indicesList=list(indicesSet)
    indicesList.sort()
    indicesSymbols=[]
    for i in indicesList:
        indicesSymbols.append(symbols(i))
    
    return res,tuple(indicesSymbols)


def indicesMultToSimp_tensor(expr):

    if isinstance(expr, Indexed):
        up_idx, down_idx = expr.indices
        return expr.base[_tupleMultTosimp(up_idx, down_idx)]
    else:
        raise TypeError("Input must be an Indexed object")

def indicesMultToSimp_mul(expr):

    res = 1
    if isinstance(expr, Mul):
        for term in expr.args:
            if isinstance(term, Indexed):
                res*= indicesMultToSimp_tensor(term)
            else:
                res*= term
    else:
        raise TypeError("Input must be a Mul object")
    
    return res

# Worker function for parallel indicesMultToSimp
# Note: _tupleMultTosimp uses global indicesSet, which may cause issues in parallel
# This worker returns (result, local_indices_set) to handle the global state
def _worker_indices_simp(term):
    """Worker function for parallel indicesMultToSimp_add."""
    global indicesSet
    indicesSet = set()  # Reset for this worker
    
    if isinstance(term, Mul):
        result = indicesMultToSimp_mul(term)
    elif isinstance(term, Indexed):
        result = indicesMultToSimp_tensor(term)
    else:
        result = term
    
    return result, indicesSet


def indicesMultToSimp_add(expr, parallel=False, show_process=True):
    """
    Simplify indices in Add expression.
    
    Args:
        expr: SymPy Add expression
        parallel: Enable parallel computing using ProcessPoolExecutor
        show_process: Show progress bar
    """
    if not isinstance(expr, Add):
        raise TypeError("Input must be an Add object")
    
    terms = list(expr.args)
    total = len(terms)
    
    if parallel and total > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        
        if show_process:
            print(f"Parallel repetitive index simplifying {total} terms...")
        
        results = []
        all_indices = set()
        
        with ProcessPoolExecutor() as executor:
            future_to_idx = {executor.submit(_worker_indices_simp, term): i 
                             for i, term in enumerate(terms)}
            
            if show_process:
                progress = ProgressBar(total, "parallel repetitive index simplifying ")
            
            for future in as_completed(future_to_idx):
                result, local_indices = future.result()
                results.append(result)
                all_indices.update(local_indices)
                
                if show_process:
                    progress.update()
        
        # Update global indicesSet with collected indices
        global indicesSet
        indicesSet = all_indices
        
        return Add(*results)
    else:
        # Sequential computing (original implementation)
        res = []
        if show_process:
            progress = ProgressBar(total, "repetitive index simplifying")
        
        for term in terms:
            if isinstance(term, Mul):
                res.append(indicesMultToSimp_mul(term))
            elif isinstance(term, Indexed):
                res.append(indicesMultToSimp_tensor(term))
            else:
                res.append(term)
            
            if show_process:
                progress.update()
        
        return Add(*res)

def indicesMultToSimp(expr, parallel=False, show_process=True):
    """
    The indices of canonicalize result are repetitive. This function is designed to simplify the indices.  
    """
    global indicesSet
    indicesSet = set()
    if isinstance(expr, Add):
        res = indicesMultToSimp_add(expr, parallel=parallel, show_process=show_process)
    elif isinstance(expr, Mul):
        res = indicesMultToSimp_mul(expr)
    elif isinstance(expr, Indexed):
        res = indicesMultToSimp_tensor(expr)
    else:
        res = expr

    indicesList=list(indicesSet)
    indicesList.sort()
    indicesSymbols=[]
    for i in indicesList:
        indicesSymbols.append(symbols(i))
    return res, tuple(indicesSymbols)

def reIndices(expr):
    """
    Re-index the expression.
    e.g. A[i,j]*B[k,l] -> A[a,b]*B[c,d]
    """
    canon_expr = qcombo.canonical.canonicalize(expr,parallel=False,show_process=False)
    reIndices_expr,indices_set =indicesMultToSimp(canon_expr,parallel=False,show_process=False)
    return reIndices_expr

################################################################
    
def get_all_indices(expr):
    """
    get all indices from expression
    """
    indices_set = set()

    for term in preorder_traversal(expr):
        if isinstance(term, Indexed):
            for i in term.indices[0]: # get up indices
                indices_set.add(i)
            for i in term.indices[1]: # get down indices
                indices_set.add(i)

    return sorted(indices_set, key=str)


def flatten_indices(nested_list):
    """
    Flatten a nested list of indices.
    """
    result = []
    for item in nested_list:
        if isinstance(item, list):
            result.extend(flatten_indices(item))
        else:
            result.append(item)
    return result

####################################################################
def find_tensor(expr, tensor = "A"):
    """
    Find the first occurrence of a tensor with the specified base name.
    
    Parameters:
        expr: SymPy expression to search
        tensor: base name string to find, default 'A'
    
    Returns:
        The first matching Indexed tensor, or None if not found
    """
    for term in preorder_traversal(expr):
        if isinstance(term, Indexed) and term.base == IndexedBase(tensor):
            return term
    return None

def antisymmetrize_tensor(expr,tensor):
    """
    Anti-symmetrize the expression with respect to the indices of a specified tensor.
    
    Parameters:
        expr: SymPy expression
        tensor: base name of the tensor to anti-symmetrize
    
    Returns:
        Anti-symmetrized expression, or original if tensor not found
    """
    To_tensor = find_tensor(expr, tensor)
    if To_tensor is not None:
        up_indices, lo_indices = To_tensor.indices[0], To_tensor.indices[1]
        return antisymmetrize(expr, up_indices, lo_indices)
    else:
        return expr

def antisymmetrize_ME(expr):
    """
    Anti-symmetrize all matrix elements (Indexed tensors) in the expression.
    
    Replaces each tensor with its anti-symmetrized form.
    
    Parameters:
        expr: SymPy expression
    
    Returns:
        Expression with all tensors anti-symmetrized
    """
    replace_dict = {}
    for term in preorder_traversal(expr):
        if isinstance(term, Indexed):
            up_indices, lo_indices = term.indices[0], term.indices[1]
            replace_dict[term] = antisymmetrize(term, up_indices, lo_indices)
    return expr.xreplace(replace_dict)          


def replace_A_tensor(expr,replace_term = 1):
    """
    Replace all occurrences of tensor A in the expression with a given term.
    
    Parameters:
        expr: SymPy expression
        replace_term: replacement value, default 1 (removes A tensor)
    
    Returns:
        Expression with all A tensors replaced
    """
    replace_dict = {}
    for term in preorder_traversal(expr):
        if isinstance(term, Indexed) and term.base == A:
            replace_dict[term] = replace_term
    return expr.xreplace(replace_dict)

def find_A_tensor(expr):
    """
    Find the first occurrence of tensor A in the expression.
    
    Parameters:
        expr: SymPy expression to search
    
    Returns:
        The first A tensor found, or None if not found
    """
    for term in preorder_traversal(expr):
        if isinstance(term, Indexed) and term.base == A:
            return term
    return None


#######################################################################
# Index permutation and anti-symmetrization functions

def swap_indices(expr, idx_a, idx_b):
    """
    Swap two indices everywhere in the expression.
    Equivalent to applying permutation operator P_{ab}.
    
    Parameters:
        expr: SymPy expression
        idx_a, idx_b: indices to swap (SymPy symbols)
    
    Returns:
        Expression with idx_a and idx_b swapped
    
    Example:
        swap_indices(A[a,b]*B[c,a], a, b)
        >> A[b,a]*B[c,b]
    """
    if idx_a == idx_b:
        return expr

    if type(idx_a) == str:
        idx_a = IndexedBase(idx_a)
    if type(idx_b) == str:
        idx_b = IndexedBase(idx_b)

    # Use a temporary symbol to avoid conflicts during swap
    tmp = symbols('__tmp_swap__')
    return expr.xreplace({idx_a: tmp, idx_b: idx_a}).xreplace({tmp: idx_b})


def antisymmetrize(expr, up_indices, lo_indices):
    """
    Anti-symmetrize the expression with respect to specified upper and lower index groups.
    
    Applies the operator: 1/(n!*m!) * prod_{pairs} (1 - P_{pair})
    
    For up_indices = [a,b] and lo_indices = [c,d]:
        result = 1/4 * (1-P_{ab})(1-P_{cd}) * expr
               = 1/4 * (expr - P_{ab}*expr - P_{cd}*expr + P_{ab}*P_{cd}*expr)
    
    This ensures the result satisfies:
        result^{ab}_{cd} = -result^{ba}_{cd} = -result^{ab}_{dc} = result^{ba}_{dc}
    
    Parameters:
        expr: SymPy expression
        up_indices: list of upper indices to anti-symmetrize, e.g., [a, b]
        lo_indices: list of lower indices to anti-symmetrize, e.g., [c, d]
    
    Returns:
        Anti-symmetrized expression (normalized by 1/(n!*m!))
    
    Example:
        expr = A[a,b,c,d]
        antisymmetrize(expr, [a,b], [c,d])
        >> 1/4 * (A[a,b,c,d] - A[b,a,c,d] - A[a,b,d,c] + A[b,a,d,c])
    """
    result = 0
    
    for up_perm in permutations(up_indices):
        for lo_perm in permutations(lo_indices):
            # Build replacement dict
            replace_dict = {}
            for old, new in zip(up_indices, up_perm):
                if old != new:
                    replace_dict[old] = new
            for old, new in zip(lo_indices, lo_perm):
                if old != new:
                    replace_dict[old] = new
            
            # Calculate sign: product of signs of upper and lower permutations
            # Permutation sign = (-1)^(number of transpositions)
            up_sign = _perm_sign(up_indices, list(up_perm))
            lo_sign = _perm_sign(lo_indices, list(lo_perm))
            sign = up_sign * lo_sign
            
            # Apply replacement and accumulate
            if replace_dict:
                term = expr.xreplace(replace_dict)
            else:
                term = expr
            
            result += sign * term
    
    # Normalize by 1/(n! * m!)
    n = len(up_indices)
    m = len(lo_indices)
    normalization = S(1) / (factorial(n) * factorial(m))
    
    return expand(result * normalization)

def antisymmetrize_expr(expr):
    """
    Anti-symmetrize the expression with respect to the indices of tensor A.
    
    Parameters:
        expr: SymPy expression containing tensor A
    
    Returns:
        Anti-symmetrized expression
    example:
        expr = (G^ai_cj H^bj_di)*A^ab_cd
        antisymmetrize_expr(expr)
        >> 1/4 * (G^ai_cd H^bj_di - G^bi_cj H^aj_di - G^ai_dj H^bj_ci + G^bi_dj H^aj_ci)*A^ab_cd
    warning: make sure only one type of A_tensor is used in the expression
    """
    A_tensor = find_A_tensor(expr)
    if A_tensor is not None:
        up_indices, lo_indices = A_tensor.indices[0], A_tensor.indices[1]
        tem_expr = expr / A_tensor
        antisymmetrized_expr = antisymmetrize(tem_expr, up_indices, lo_indices)
        return antisymmetrized_expr * A_tensor
    else:
        return expr

def _perm_sign(original, permuted):
    """
    Calculate the sign of a permutation.
    
    Parameters:
        original: list in original order
        permuted: list in permuted order
    
    Returns:
        +1 for even permutation, -1 for odd permutation
    """
    n = len(original)
    if n <= 1:
        return 1
    
    # Map to integer indices
    mapping = {v: i for i, v in enumerate(original)}
    perm_list = [mapping[v] for v in permuted]
    
    # Count inversions
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if perm_list[i] > perm_list[j]:
                inversions += 1
    
    return 1 if inversions % 2 == 0 else -1


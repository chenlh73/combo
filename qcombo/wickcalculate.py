'''
This module includes operations for handling Wick's theorem cases.
'''
from sympy.utilities.iterables import multiset_partitions  # For set partitions
from itertools import permutations  # For permutations
from sympy import IndexedBase  # For storing indices
from qcombo.tools import ProgressBar
from sympy.core.add import Add



# Alpha=IndexedBase(chr(913))
Alpha=IndexedBase('A')
xi=IndexedBase(chr(958))
lamda=IndexedBase(chr(955))
delta=IndexedBase(chr(948))
n=IndexedBase('n')

class _ListOperations:
    '''
    Some important operations to process the list. 
    
    '''

    def transpose(listObject)->list:
        '''
        Transpose a list.

        '''
        #Use numpy to trans it may be more convenient.But in This way, we dont need extra package.
        res=list(map(list,zip(*listObject)))
        return  res
    
    def signature(listObject)->int:
        '''
        Make sure that the listObjetct is a 1 dimension list.
        If the list is odd permutation,return -1;Otherwise ,return 1.
        '''
        sig=1
        for i in range(len(listObject)):
            for listObject_j in listObject[:i]:
                if(str(listObject_j)>str(listObject[i])):
                    sig=(-1)*sig
        return sig
    
    def sort(listObject)->list:
        '''
        Sort the list by the length of it element. If the element is list, sort the listObject in the length of it's sublist. 
        '''
        listObject.sort(key=len,reverse=False)
        #lambda listObject:listObject.sort(key=len,reverse=False)
        return listObject

    def union(listObject)->list:
        '''
        remove the same elements of listObject and keep one of them in the listObject.
        '''
        res=[]
        for listObject_i in listObject:
            if listObject_i not in res:
                res.append(listObject_i)
        return res

    def setPartitions(listObject)->list:
        '''
        Return unique partitions of the given multiset (in list form)

        # Examples
        >>> setPartitions([1, 1, 2])
        [[[1, 1, 2]], [[1, 1], [2]], [[1, 2], [1]], [[1], [1], [2]]]

        # Note
        Used `multiset_partitions` in  `sympy.utilities.iterables` to achieve setPartitions. 

        '''
        res=[]
        for i in multiset_partitions(listObject):# trans multiset_partitions object to list
            res.append(i)
        return res

    def sortPartitions(listObject)->list:
        '''
        Sort the list by the length of it's elements.Meanwhile sort it elemnts.

        # Examples
        >>> sortPartitions([[[1,2],[1,2,3],[1]],[[1,2,3,4],[1,2,3]]])
        [[[1, 2, 3], [1, 2, 3, 4]], [[1], [1, 2], [1, 2, 3]]]
        >>> sortPartitions([[[1,2],[1,2,3],[1]]])
        [[[1], [1, 2], [1, 2, 3]]]    

        '''
        res=listObject
        if(len(res)>=1):
            res.sort(key=len,reverse=False)
            for i in res:
                i.sort(key=len,reverse=False)
        return res 
 

# Worker function for parallel computing (must be at module level for pickling)
def _worker_construct_terms(args):
    """Worker function for ProcessPoolExecutor in generalizedWick."""
    part_i, idxUp, idxLo, LUp, LLo, RUp, RLo, wickMode = args
    wick = Wick([LUp, LLo], [RUp, RLo])
    wick.wickMode = wickMode  # Inherit parent process's wickMode
    return wick._ConsturctTerms(part_i, idxUp, idxLo, LUp, LLo, RUp, RLo)


class Wick:
    '''
    This class is used to calculate the Wick's theorem.
    '''
    def __init__(self,LUpLo,RUpLo):
        self._LUp,self._LLo,self._RUp,self._RLo=LUpLo[0],LUpLo[1],RUpLo[0],RUpLo[1]
        self.gw=0#generalizedWick result
        self.cmt=0#commutate result

        # clh, add wickMode attribute to control the wick theorem mode, 'SR' or 'MR'
        self.wickMode = 'MR'# 'SR' or 'MR',single-reference or multi-reference
        self.wickShowProcess = True# Display the process bar of calculation
        self.wickParallel = False  # Enable parallel computing using ProcessPoolExecutor
    
    def _MatchPartitions(self,pUp,pLo):
        '''
        Match partitions where the upper and lower index lengths are equal.
        The sub-list lengths must also match to form valid index pairings.
        '''
        res=[]
        for pUp_i in pUp:
            lengthPartitionUp=len(pUp_i)
            for pLo_j in pLo:
                lengthPartitionLo=len(pLo_j)
                if(lengthPartitionUp==lengthPartitionLo):
                    strcutPartitionUp=list(map(len,pUp_i))#return lengths of sub-lists
                    strcutPartitionLo=list(map(len,pLo_j))
                    if(strcutPartitionUp==strcutPartitionLo):#sub-list lengths must also be identical
                        res.append(_ListOperations.transpose([pUp_i,pLo_j]))   
        return res


    # def _EvaluateConstraction(self,full_i_k,LUp,LLo,RUp,RLo):
    #     '''
    #     evaluate a contraction for a specific combination of indices
    #     '''
    #     if(((set(full_i_k[0]) & set(LUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(LLo)==set(full_i_k[1])))  or  ((set(full_i_k[0]) & set(RUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(RLo)==set(full_i_k[1])))):
    #         res=0
    #     else:
    #         if(len(full_i_k[0])==1 and len(full_i_k[1])==1  and  (set(full_i_k[0]) & set(RUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(LLo)==set(full_i_k[1]))):
    #             #res=symbols_laxexp('xi',full_i_k))
    #             res=xi[tuple(full_i_k[0]),tuple(full_i_k[1])]
    #         else:
    #             #res=symbols_laxexp('lambda',full_i_k) 
    #             res=lamda[tuple(full_i_k[0]),tuple(full_i_k[1])]
    #     return res

    def _EvaluateConstraction(self,full_i_k,LUp,LLo,RUp,RLo):
        '''
        evaluate a contraction for a specific combination of indices
        '''
        #modified by clh
        # 
        if self.wickMode == 'MR': # MR mode: multi-body lambda is allowed
            if(((set(full_i_k[0]) & set(LUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(LLo)==set(full_i_k[1])))  or  ((set(full_i_k[0]) & set(RUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(RLo)==set(full_i_k[1])))):
                res=0
            else:
                if(len(full_i_k[0])==1 and len(full_i_k[1])==1  and  (set(full_i_k[0]) & set(RUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(LLo)==set(full_i_k[1]))):
                    #res=symbols_laxexp('xi',full_i_k))
                    res=xi[tuple(full_i_k[0]),tuple(full_i_k[1])]
                else:
                    #res=symbols_laxexp('lambda',full_i_k) 
                    res=lamda[tuple(full_i_k[0]),tuple(full_i_k[1])]
        
        elif self.wickMode == 'SR': # SR mode: only single-body lambda is allowed
            if(((set(full_i_k[0]) & set(LUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(LLo)==set(full_i_k[1])))  or  ((set(full_i_k[0]) & set(RUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(RLo)==set(full_i_k[1])))):
                res=0 # All contraction indices come from a single operator, cannot contract
            else:
                if(len(full_i_k[0])==1 and len(full_i_k[1])==1  and  (set(full_i_k[0]) & set(RUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(LLo)==set(full_i_k[1]))):
                    res=xi[tuple(full_i_k[0]),tuple(full_i_k[1])] # xi: contraction between upper-right and lower-left indices

                elif (len(full_i_k[0])==1 and len(full_i_k[1])==1  and  (set(full_i_k[0]) & set(LUp)==set(full_i_k[0])) and (set(full_i_k[1]) & set(RLo)==set(full_i_k[1]))):
                    res=lamda[tuple(full_i_k[0]),tuple(full_i_k[1])] #lambda: contraction between upper-left and lower-right indices

                else:
                    res=0 # All other cases

        return res


    def _ConsturctTerms(self,part_i,idxUp,idxLo,LUp,LLo,RUp,RLo):#part_i: one specific upper/lower index combination, i.e. one output configuration
        '''
        '''
        # e.g. part[i] = [[[a], [b]], [[p], [r]]]
        sigUp=_ListOperations.signature(idxUp)#determine original sign ordering, since input indices may not be in ascending order
        sigLo=_ListOperations.signature(idxLo)
        #Take lower indices after transpose, and enumerate all permutations of lower indices
        #Reason: combinations from Step are incomplete because indices were sorted ascendingly,
        #but actual combinations may not follow ascending order.
        permutations_part_i=list(map(list,list(permutations(_ListOperations.transpose(part_i)[1]))))#using itertools.permutations for full permutation
        #Re-sort and deduplicate: only subsets of same length are permuted
        permutations_part_i=_ListOperations.union(map(_ListOperations.sort,permutations_part_i))#sort and deduplicate
        full=list(map(#convert map output to list
            _ListOperations.transpose,list(#transpose again, one sub-list per upper/lower group#this list may be unnecessary
            map(#combine each lower-index permutation with upper indices
            lambda permu_i:[_ListOperations.transpose(part_i)[0],permu_i],#transpose part_i to get upper indices, then combine with lower
            permutations_part_i)
            )
        ))

        # full contains all possible permutations/ combinations for the same part_i
        # e.g. part[i] = [[[a], [b]], [[p], [r]]]
        # full =  [[[[a], [b]], [[p], [r]]], [[[a], [r]], [[p], [b]]]]
        expr=0
        for i in range(len(full)):
            termIndexUp= [y for x in [y[0] for y in full[i]] for y in x]#upper index set of i-th combination: take all upper indices from full[i] then flatten
            #termIndexUp=list(np.array(full[i],dtype=object)[:,0].flatten())
            ##This would work if full were a numpy array, but avoided to reduce dependencies#This line gives wrong result because flatten cannot handle arrays with unequal sub-list lengths
            termIndexLo= [y for x in [y[1] for y in full[i]] for y in x]
            
            tmplambda=1
            for j in range(len(full[i])):#len(full[i]) = number of operator groups?
                tmpA=1
                for k in range(len(full[i])):
                    if(j==k):
                        #tmpA*=symbols_laxexp('A',full[i][k])#map indices to corresponding A operator
                        tmpA*=Alpha[tuple(full[i][k][0]),tuple(full[i][k][1])]
                    else:
                        tmpA*=self._EvaluateConstraction(full[i][k],LUp,LLo,RUp,RLo)#for other indices: xi and lambda
                expr+=_ListOperations.signature(termIndexUp)*_ListOperations.signature(termIndexLo)*tmpA
                tmplambda*=self._EvaluateConstraction(full[i][j],LUp,LLo,RUp,RLo)
            expr+=_ListOperations.signature(termIndexUp)*_ListOperations.signature(termIndexLo)*tmplambda

        expr*=sigUp*sigLo
        return expr  

    
    def generalizedWick(self):
        """
        Perform generalized Wick theorem expansion.
        
        Enumerates all valid index partitions and contractions,
        generating terms with A operators, xi (contraction) and
        lambda (density matrix) tensors.
        
        Result is stored in self.gw.
        """
        idxUp=self._LUp + self._RUp
        idxLo=self._LLo + self._RLo
        partUp=_ListOperations.sortPartitions(_ListOperations.setPartitions(idxUp))
        partLo=_ListOperations.sortPartitions(_ListOperations.setPartitions(idxLo))
        part=self._MatchPartitions(partUp,partLo)#part: all valid upper/lower index partition pairs with matching subset lengths

        if self.wickParallel:
            # Parallel computing using ProcessPoolExecutor
            from concurrent.futures import ProcessPoolExecutor, as_completed
            
            # Prepare task arguments
            task_args = [
                (part_i, idxUp, idxLo, self._LUp, self._LLo, self._RUp, self._RLo, self.wickMode)
                for part_i in part
            ]
            
            total_tasks = len(part)
            if self.wickShowProcess:
                print(f"Parallel computing with {total_tasks} tasks...")
            
            # Execute in parallel with progress display
            results = []
            with ProcessPoolExecutor() as executor:
                # Submit all tasks
                future_to_idx = {executor.submit(_worker_construct_terms, args): i 
                                 for i, args in enumerate(task_args)}
                
                # Collect results with progress bar
                if self.wickShowProcess:
                    progress = ProgressBar(total_tasks, "parallel wick")
                
                for future in as_completed(future_to_idx):
                    results.append(future.result())
                    if self.wickShowProcess:
                        progress.update()
            
            # Sum up all results
            self.gw += Add(*results)
            
            if self.wickShowProcess:
                print("Parallel computing completed!")
        else:
            # Sequential computing (original implementation)
            if self.wickShowProcess: progress = ProgressBar(len(part), "generalize wick caculating") # Display the process bar of calculation
            
            # e.g. for operator A^a_b * B^p_r, part = [[[[a, p], [b, r]]], [[[a], [b]], [[p], [r]]]]
            for i in range(len(part)):
                # e.g. part[i] = [[[a], [b]], [[p], [r]]]
                self.gw += self._ConsturctTerms(part[i],idxUp,idxLo,self._LUp,self._LLo,self._RUp,self._RLo)
                
                if self.wickShowProcess: progress.update() # Display the process bar of calculation


    def commutator(self):
        '''
        Calculate the commutator
        '''
        LUp,LLo,RUp,RLo=self._LUp,self._LLo,self._RUp,self._RLo

        positiveWick=Wick([LUp,LLo],[RUp,RLo])
        positiveWick.wickMode = self.wickMode 
        positiveWick.wickShowProcess = self.wickShowProcess
        positiveWick.wickParallel = self.wickParallel
        positiveWick.generalizedWick()
        
        negativeWick=Wick([RUp,RLo],[LUp,LLo])
        negativeWick.wickMode = self.wickMode
        negativeWick.wickShowProcess = self.wickShowProcess
        negativeWick.wickParallel = self.wickParallel
        negativeWick.generalizedWick()

        self.gw=positiveWick.gw
        self.cmt=positiveWick.gw-negativeWick.gw




if __name__=='__main__':
    # test

    a=IndexedBase('a')
    b=IndexedBase('b')
    c=IndexedBase('c')
    d=IndexedBase('d')
    e=IndexedBase('e')
    f=IndexedBase('f')
    p=IndexedBase('p')    
    q=IndexedBase('q')
    s=IndexedBase('s')
    r=IndexedBase('r')
    u=IndexedBase('u')
    v=IndexedBase('v')
    from sympy import false, true

    test=Wick([(a,b),(d,e)],[(p,q,u),(r,s,v)])
    # test=Wick([(a,b,c),(d,e,f)],[(p,q,u),(r,s,v)])
    # test=Wick([(a,),(d,)],[(p,),(r,)])

    my_wick_test = test

    my_wick_test.wickParallel = True
    my_wick_test.commutator()

    parallel_result = my_wick_test.cmt
    print("Parallel result finished")

    my_wick_test=test
    my_wick_test.wickParallel = False
    my_wick_test.commutator()

    no_parallel_result = my_wick_test.cmt
    print("No-parallel result finished")

    print(parallel_result == no_parallel_result)

    res=my_wick_test.cmt
    # print(res)


#TODO:
#1.wick_mode = "SR" should be improved to decrease time-consuming
# now the "SR" seems using the same time as "MR"

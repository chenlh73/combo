if __name__ == "__main__":
    import qcombo

    left = 2
    right = 3
    contraction = None # None means all contraction
    latexOutput = "./examples/results/commutator_{}B{}B.tex".format(left, right)
    amcOutput = "./examples/results/commutator_{}B{}B.amc".format(left, right)
    wick_mode  = "MR" # 'SR' or 'MR', single-reference or multi-reference
    # for larger computations like [2,3], enable parallel computing to save time
    parallel = False # Enable parallel computing using ProcessPoolExecutor
    show_process = True # Display the process bar of calculation
    savefile = True # Save the result to file

    qcombo.easyCombo(left=left,right=right,latexOutput=latexOutput,
                        amcOutput=amcOutput, wick_mode=wick_mode, parallel=parallel, 
                        show_process=show_process, savefile=savefile)




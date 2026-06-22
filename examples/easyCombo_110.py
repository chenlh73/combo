if __name__ == "__main__":
    import qcombo

    wick_mode  = "MR" # 'SR' or 'MR', single-reference or multi-reference
    parallel = False # Enable parallel computing using ProcessPoolExecutor
    show_process = True # Display the process bar of calculation
    savefile = True # Save the result to file

    commutator_110 = qcombo.easyCombo(left=1,right=1,contraction=0,latexOutput="./examples/results/test_commutator_110.tex",
                        amcOutput="./examples/results/test_commutator_110.amc", wick_mode=wick_mode, parallel=parallel, 
                        show_process=show_process, savefile=savefile)

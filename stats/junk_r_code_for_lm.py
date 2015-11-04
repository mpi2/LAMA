__author__ = 'neil'
  tstats = []
        pvals = []

        # np.array_split provides a split view on the array so does not increase memory
        # The result will be a bunch of arrays split down the second dimension


        # write csv to tempfile for R
        raw_data_file = join(tempfile.gettempdir(), 'raw_data_for_r.csv')
        groups_file = join(tempfile.gettempdir(), 'groups_for_liear_model.csv')

        groups = ['wildtype'] * len(self.wt_data)
        groups.extend(['mutant'] * len(self.mut_data))
        group_name = 'sex'

        r_pvals_out = join(tempfile.gettempdir(), 'r_pvals_out.csv')
        r_tvals_out = join(tempfile.gettempdir(), 'r_tvals_out.csv')

        # Save the group info.
        with open(groups_file, 'wb') as gf:
            gf.write(group_name + '\n')
            for group in groups:
                gf.write(group + '\n')

        chunked_mut = np.array_split(self.mut_data, 300, axis=1)
        chunked_wt = np.array_split(self.wt_data, 300, axis=1)

        script_dir = os.path.dirname(os.path.realpath(__file__))
        r_linear_model_script = join(script_dir, LINEAR_MODEL_SCIPT)

        pval_results = []
        tval_results =[]

        for wt_chunks, mut_chunks in zip(chunked_wt, chunked_mut):
                # Write the chunk of raw data

            with open(raw_data_file, 'ab') as rf:
                np.savetxt(rf, wt_chunks, delimiter=',', fmt='%10.3f')  # Maybe increase decimal places
                np.savetxt(rf, mut_chunks, delimiter=',', fmt='%10.3f')

            subprocess.check_output(['Rscript',
                                     r_linear_model_script,
                                     raw_data_file,
                                     groups_file,
                                     r_pvals_out,
                                     r_tvals_out])

            # read in the results
            pval_results.append(np.genfromtxt(r_pvals_out, delimiter=','))
            tval_results.append(np.genfromtxt(r_tvals_out, delimiter=','))

            try:
                os.remove(r_pvals_out)
                os.remove(r_tvals_out)
                os.remove(raw_data_file)
            except OSError:
                pass

        pvals = np.hstack(pval_results)
        tvals = np.hstack(tval_results)
        print pvals
        print tvals

                # make sure to delete temp file after doing lm on a chunk

            # Run R script to do linear model analysis

            #subprocess.check_output(['Rscript', ])


        tstats_chunk, pval_chunk = self.runttest(wt_chunks, mut_chunks)
        pval_chunk[np.isnan(pval_chunk)] = 0.1
        pval_chunk = pval_chunk.astype(np.float32)
        tstats.extend(tstats_chunk)
        pvals.extend(pval_chunk)

        pvals = np.array(pvals)
        tstats = np.array(tstats)


        fdr = self.fdr_class(pvals)
        qvalues = fdr.get_qvalues(self.mask)
        gc.collect()

        #self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues)
        self.filtered_tscores = self._result_cutoff_filter(tstats, qvalues) # modifies tsats in-place

        # Remove infinite values
        self.filtered_tscores[self.filtered_tscores > MINMAX_TSCORE] = MINMAX_TSCORE
        self.filtered_tscores[self.filtered_tscores < -MINMAX_TSCORE] = - MINMAX_TSCORE

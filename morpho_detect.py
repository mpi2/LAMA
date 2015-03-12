














if __name__ == '__main__':

    import argparse


    parser = argparse.ArgumentParser("messageS")
    parser.add_argument('-c', dest='cubesize', type=int, help='Size of the sub array', required=True)
    parser.add_argument('-d', dest='deforms', help='Folder containing deformation field files ', default=None)
    parser.add_argument('-j', dest='jacobians', help='Folder containing jacobian field files ', default=None)
    parser.add_argument('-r', dest='registered_vols', help='Folder containing registered vols, for intensity difference'
                                                           ' analysis', default=None)
    parser.add_argument('-o', dest='outfile', help='File to store pickle file of means/stdvs of vectors,intensities etc', required=True)
    parser.add_argument('-t', dest='threads', type=int, help='How many threads to use', default=4)
    parser.add_argument('-dp', dest='deform_file_pattern', help='String that is contained in the deform file names',
                        default='deformationField')
    parser.add_argument('-jp', dest='jac_file_pattern', help='String that is contained in the jacobian file names',
                        default='spatialJacobian')

    args = parser.parse_args()
    if not args.deforms and not args.jacobians and not args.registered_vols:
        sys.exit('You need to supply deformations, jacobians, or registered volumes (-d, -j, -r)')

    #cProfile.run(run(args.deforms, args.cubesize, args.outfile))
    global chunksize
    chunksize = args.cubesize
    run(args)


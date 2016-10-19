infile <- '/tmp/tempdata.dat'
con <- file(infile, "rb")
dim <- readBin(con, "integer", 2)
Mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
close(con)

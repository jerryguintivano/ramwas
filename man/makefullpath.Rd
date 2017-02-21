\name{makefullpath}
\alias{makefullpath}
\title{
    Combine Path and Filename into  Filename with Path
}
\description{
    Combine a path with a filename into filename with path.
}
\usage{
makefullpath(path, filename)
}
\arguments{
    \item{path}{
        Path, relative to which the filename is expected to be.
        Can be absolute, relative, or \code{NULL}.
}
    \item{filename}{
        Can be just filename, include relative path, or
        include absolute path.
}
}
\details{
    Function returns \code{filename} if it includes absolute path or
    if \code{path} is \code{NULL}.
}
\value{
    Filename with the path included.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See also \code{\link{isAbsolutePath}}.
}
\examples{
makefullpath("dir1/dir2", "file.txt")
# "dir1/dir2/file.txt"

makefullpath("dir1/dir2", "dir3/file.txt")
# "dir1/dir2/dir3/file.txt"

# Path is ignored if the filename already includes absolute path

makefullpath("dir1/dir2", "/file.txt")
# "/file.txt"

makefullpath("dir1/dir2", "C:/file.txt")
# "C:/file.txt"
}

###########################################################################
# Statistics for Microarray Analysis
# Misc. functions
#
# Date : March 19, 2001
#
# Authors: Yee Hwa (Jean) Yang and Jessica Mar
##########################################################################

########################################################################/** 
# \name{is.odd}
# \alias{is.odd}
# \alias{is.even}
# 
# \title{ Determining if a Value is Odd or Even }
# \description{
# A logical flag which determines if a value supplied by the user is
# odd or even. 
# }
# \usage{
# is.odd(x)
# is.even(x)
# }
# 
# \arguments{
#  \item{x}{integer value}
# }
# 
# }
# \value{Logical vectors \code{TRUE} or \code{FALSE} are returned
#   depending on whether the value is odd or even.
# 
# }
# 
# \author{ Jessica Mar }
# 
# \examples{
# is.odd(4)
# ## FALSE
# is.even(100)
# ## TRUE
# }
# 
# \keyword{odd, even}
#*/########################################################################

is.even <- function(x)
  {if(is.numeric(x))
     {if (x %% 2==0) {TRUE}
      else {FALSE}
     }
   
   else{
       print("Warning: Input must be an integer value")
     } 
 }

is.odd <- function(x)
  {if(is.numeric(x))
    {if (x %% 2 == 0) {FALSE}
     else {TRUE}}
 
   else{
       print("Warning: Input must be an integer value")
     }
 }

########################################################################/** 
# \name{id2image}
# 
# \alias{id2image}
# \alias{image2id}
#  
# \title{Converting an id tag to a Set of Image Coordinates and Vice Versa}
# 
# \description{
#  The function \code{id2image} converts an id tag of a gene supplied
#  by a user into a set of image coordinates regarding the location of
#  the gene being 
#  considered. Conversion of image coordinates to an id tag is performed
#  by \code{image2id}.
# }
# 
# \usage{
# id2image(X, layout)
# image2id(x, layout)
# }
# 
# \arguments{
#  \item{X}{an integer value representing the id of a particular gene}
#  \item{layout}{}
#  \item{x}{a vector of 4 integer elements which make up the image
#    coordinates of the gene.}
# }
# \details{
# The image coordinates of a gene correspond to the gene's grid row and
# grid column position within a slide, and the gene's row and column
# position within a grid.  
# }
# 
# \value{
#   \code{id2image} returns a vector of 4 integer elements which is
#   the set of image coordinates. 
#   \code{image2id} returns an integer value which is the gene's id tag.
# }
# 
# \author{Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu}}
# 
# \seealso{\code{\link{MouseArray}}}
# 
# \examples{data(MouseArray)
# # mouse.setup <- init.grid()
# 
# id2image(1024, mouse.setup)
# ## You will see: [1]  1 3 11 16
# ## the grid in which gene 1024 can be found, is in row 1, column 3
# ## and the gene is located in row 11, column 16 of this particular grid.
#  
# image2id(c(2,4,6,8), mouse.setup)
# ## You will see: [1] 2906
# ## the gene located in row 6, column 8 in the grid that is in row 2 and
# ## column 4 is the 2906th gene of the data set.  
# }
# 
# \keyword{image, id}
#*/#######################################################################
 
id2image <- 
function(X, layout)
{
        Grid.row <- layout$ngrid.r; Grid.col <- layout$ngrid.c
        Spot.row <- layout$nspot.r; Spot.col <- layout$nspot.c

        Coord <- rep(0, 4)
        Spot.size <- Spot.row * Spot.col
        ## Calculate Grid row & column coordinates
        Coord[1] <- ((X - 1) %/% (Grid.col * Spot.size)) + 1
        count.Spots <- ((X - 1) %/% Spot.size) + 1
        Coord[2] <- ((count.Spots - 1) %% Grid.col) + 1
        ## Calculate Spot row & column coordinates
        Spot.pt <- X - (count.Spots - 1) * Spot.size
        Coord[3] <- ((Spot.pt - 1) %/% Spot.col) + 1
        Coord[4] <- ((Spot.pt - 1) %% Spot.col) + 1
        Coord
}

image2id <- 
function(x, layout)
{
        Grid.row <- layout$ngrid.r; Grid.col <- layout$ngrid.c
        Spot.row <- layout$nspot.r; Spot.col <- layout$nspot.c

        temp <- Spot.col * Spot.row
        temp * ((x[1] - 1)*Grid.col+(x[2] - 1))+(x[3] - 1)*Spot.col+ x[4]
}



##########################################################################
#                                End of file
##########################################################################

estimateLAI <- function(THETA,GAP,option=list(),...)
{
  if (!hasArg("THETA") || !hasArg("GAP"))
  {
    warning("Something is wrong with the input. Make sure the first and second inputs are provided!")
    return(NULL)
  }

  if ( !(is.vector(THETA) || is.matrix(THETA) ))
  {
    warning("Something is wrong with the inputs. Make sure THETA and GAP are either vectors or matrices")
    return(NULL)
  }

  if ( is.vector(THETA) && !( is.vector(GAP) && length(THETA)==length(GAP) ) )
  {
    warning("Something is wrong with the inputs. Make sure THETA and GAP are vectors of the same length")
    return(NULL)
  }

  if ( is.matrix(THETA) && !( is.matrix(GAP) && dim(THETA)[1]==dim(GAP)[1] && dim(THETA)[2]==dim(GAP)[2])  )
  {
    warning("Something is wrong with the inputs. Make sure THETA and GAP are matrices of the same size")
    return(NULL)
  }
  
  if (!is.list(option))
  {
    warning("Something is wrong with the third parameter. Make sure it is a list.")
    return(NULL)
  }

 
 
  
  #if(season)
  #{
  #res=.Call(SARAH_beastST_multipleChain_fast, data, OPTION)
  #}
  #else
  #{
  #res=.Call(SARAH_beastTrend_multipleChain_fast, data, OPTION)
  #}
  

  ANS=.Call("LAI_estimate", THETA,GAP,option)
  #cat(names(ANS))
  class(ANS)=c("hemiphoto","data.frame")
  return(ANS)
  #cat("\n\n==================================================================\n")
  #s=paste("The beast output variable (e.g., x) is a LIST object. Type names(x)", 
  #"to see a list of elements in x.  \n\nThe current 'x' contains the following #elements: ")
  #cat(s)
  #namesList=names(ANS);
  #cat(namesList)
  #cat(". ")
  #cat("\n\nCheck individual elements to see the model outputs (e.g,type x$t to see #the fitted trend) or plot(x) to draw the model decomposition result.")
 #cat("\n==================================================================")
	
  
 
 
  
    
}
 

 
 
 

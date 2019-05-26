simHemiphoto <- function(ladType=1,ladParameter, option=list(),...)
{
  if (!hasArg("ladType") || is.list(ladType))
  {
    warning("Something is wrong with the input. Make sure the first input is a string or integer to specify the LAD")
    return(NULL)
  }

if (!hasArg("ladParameter") || length(ladParameter)==1)
  {
    warning("Something is wrong with the second parameter. Make sure the second paramer  is a vector of two numbers.")
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
  

  ANS=.Call("LAI_simulate", ladType,ladParameter, option)
  cat(names(ANS))
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
 

 
 
 

print.hemiphoto<-function (x,...)
{
#https://stackoverflow.com/questions/50561768/r-get-argument-names-from-function-call
  #argName=formalArgs(deparse(substitute(x))[[1]])
#https://stackoverflow.com/questions/4959463/getting-the-object-name-for-s3-print-method-failing
 #argName=deparse(substitute(x))
 
 #argName=strtrim(argName)
  #nTS= length(x$tN);
  print.data.frame(x, digits=4);
   
  cat("\n\n==================================================================\n")
  s=paste("As you can see from the above, the estimateLAI output variable (e.g., x) is a data frame object. Type names(x)", 
  "to see a list of elements in x.  \n\nThe current 'x' contains the following elements: ")
  cat(s)
  namesList=names(x);
  cat(namesList)
  cat(". ")
  cat("\n\nCheck individual elements to see the model outputs (e.g, type x$LAI_BNR to see the estimated LAI values based on 19 different LAD models).")
   cat("\n==================================================================")
   
   
  
}

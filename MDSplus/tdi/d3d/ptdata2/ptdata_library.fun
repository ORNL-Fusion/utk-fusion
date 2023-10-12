/**************************************************************************************************************************
        PUBLIC FUN PTDATA_LIBRARY()
        
        This TDI function retrieves the location of the PTDATA Library from the PTDATA_LIBRARY environment variable

        Author:         Sean Flanagan (flanagan@fusion.gat.com) 

**************************************************************************************************************************/

PUBLIC FUN PTDATA_LIBRARY() {
	_lib = TranslateLogical("PTDATA_LIBRARY");
	RETURN(_lib);
}

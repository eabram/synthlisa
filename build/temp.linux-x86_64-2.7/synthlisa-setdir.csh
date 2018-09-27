if !($?PYTHONPATH) then
    setenv PYTHONPATH 
else
    setenv PYTHONPATH :$PYTHONPATH
endif


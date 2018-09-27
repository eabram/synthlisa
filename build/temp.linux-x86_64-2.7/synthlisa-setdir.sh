if [ -z "${PYTHONPATH}" ]
then
    PYTHONPATH=""; export PYTHONPATH
else
    PYTHONPATH=":$PYTHONPATH"; export PYTHONPATH
fi


#htdecodetoken
#ID=htdecodetoken | awk -F '/' '{print $4}'
#ID="51596"
#echo "id -u = ${ID}"
#htgettoken -i minerva -a htvaultprod.fnal.gov:8200
#htdecodetoken
#USER_PROXY="/tmp/x509up_minerva_Analysis_${ID}"
#export X509_USER_PROXY=$USER_PROXY
#echo $X509_USER_PROXY
#setup jobsub_client v_lite
#voms-proxy-destroy;
#kx509;
#voms-proxy-init -rfc --voms=fermilab:/fermilab/minerva/Role=Analysis --noregen -valid 24:0;

export EXPERIMENT=minerva
export IFDH_DEBUG=0

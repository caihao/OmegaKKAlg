# echo "Setting OmegaKKAlg OmegaKKAlg-00-00-01 in /panfs/panfs.ihep.ac.cn/home/data/caih/6.5.5/mywork"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/tool/CMT/v1r20
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=OmegaKKAlg -version=OmegaKKAlg-00-00-01 -path=/panfs/panfs.ihep.ac.cn/home/data/caih/6.5.5/mywork  -no_cleanup $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}


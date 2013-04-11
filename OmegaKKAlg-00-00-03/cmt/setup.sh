# echo "Setting OmegaKKAlg OmegaKKAlg-00-00-01 in /panfs/panfs.ihep.ac.cn/home/data/caih/6.5.5/mywork"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/ihep.ac.cn/bes3/offline/tool/CMT/v1r20; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh

tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=OmegaKKAlg -version=OmegaKKAlg-00-00-01 -path=/panfs/panfs.ihep.ac.cn/home/data/caih/6.5.5/mywork  -no_cleanup $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}


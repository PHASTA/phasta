c..............................................................................

        subroutine findIsrfid(isrfid)
          
        include "common.h"

        integer globNodemap(nshg)
        integer nnodesFound
        integer isrfid(nshg)
        isrfid=-1

        do isfID =0, MAXSURF
          call sfID2np(isfID,nnodesFound,globNodemap)  
          if(nnodesFound .gt. 0)then
           do i=1,nnodesFound
             nn=blobNodemap(i)
             if(isrfid(nn).eq.-1) then
                 isrfid(nn)=isfID
             elseif(isrfid(nn).ne.isfID) then ! different srfid on shared model edge or vertex should get set to a value we know to ignore.  This allowsit to be KEPT if same 
                 isrfid(nn)=-2
             endif
           enddo
          endif
        enddo
         
        return 
        end

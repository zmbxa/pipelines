name=$1
date=$2

show_usage="Usage: bash download.obs.sh <name> <date> \n"

if [[ -z ${name} || -z ${date} ]]; then
        echo -e $show_usage
        exit 0
fi

case "$name" in 
                HSQ) I=U59XDM5IPA1LHM3XW58X; K=fAD84AlCSf1nutTuaO0jjef0hd7O3F29zQANNPq3; O=westlakehousq;;
                MWY) I=KMA6BAH989S0ZDR4SSHS; K=pBwVCx2bB5yTqthQCBcWdRAYSaueuz7uWeUWB8b3; O=westlakemwy;;
                GY)  I=L2JB6YPDDCSYYOT5ADDT; K=ew7Ue0x1OH4QfDZYiy3oEv3nuCQiL4i3umkUc4Vw; O=westlakeguoyi;;
                HJC) I=CBN8ND08KRC1CNDTQXVI; K=AQe7j2benUltCEb8ZcNQPa8ZrphYeDEg2YRbro5H; O=westlakehoujc;;
                HJH) I=0B9UYBALXW1CCNEGJULV; K=TPpFXOWKYFYTUnqFWZSG23ai9G8wNBGx0CuiLojk; O=fudanyulab;;
                LLX) I=P5XWXOMBKIDUNWIJIKM7; K=bk89uGM4Huz5E29T1EPUiqbgL03ZAQVRaPOhlnG9; O=westlakelilx;;
                WX)  I=POIRBO8UKIB8FNDNTAKJ; K=wFFmxheppibeNJ9D7DrVAOU5ueFN04lgktlsGtlP; O=westlakewangxiao;;
                XX)  I=JF8HMBQFYAA3L6AY71DE; K=USEa5OqiARa9vMy66lZ7qmRwCM2QR755abtpmlR6; O=westlakexiongx;;
                HM)  I=2LSQLCMCZ5265J7MBBIJ; K=jH8aIRhNcvg2UcmwuvxUNvLh4brdb5m1RbjaEGYq; O=westlakehuangm;;
                CYQ) I=LA1KF1ONRNZMKHHK6952; K=9DlhE4WtiwbycUq0HLau7PzcJzY4qKclnQJYqECf; O=westlakechenyq;;
                NTY) I=BULQM7BVM9B2PIJUS7FD; K=XxVQD2MfbKNDT24tOTOOEtvzvXMgXb0aJlx1mxBh; O=westlakenity;;
                SYJ) I=81FTB0UWWBF8RWTB8UHK; K=jTMuev2YcGbq42lV0JPWE4aCFgklNnHeZY7vOrIx; O=westlakeshiyj;;
                SZJ) I=ZAZO6E6HAYQLY5AYBCJW; K=mEI70asgiQQ80uaNtjzybFMbBbbZaXeIWuPt3vvr; O=westlakesuzj;;
                ZYX) I=XL0LY70ZRY3AXMSPL1KG; K=asnttF8bduakqofJJzfDAz1WJiMDBbmS0gVUSZaA; O=westlakezhangyx;;
                NYX) I=G2AGPNFBFS2WEAVBAODT; K=d0Vi8oGU50le1xOjWn8IkUIRuIbuy0EkgS9H198H; O=westlakeniuyuxiao;;
                *)
                esac

obsutil config -i=$I -k=$K -e=http://218.94.125.245:8999

obsutil cp obs://$O/$date ./$name -f -r -u


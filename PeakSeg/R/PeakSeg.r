PeakSeg<-function(vec,winnum=3,step=20){
quan = function(oo, alpha){
	a1 = as.numeric(quantile(sort(oo), probs = c(alpha, 1 - alpha))[1])
	b1 = as.numeric(quantile(sort(oo), probs = c(alpha, 1 - alpha))[2])
	result = oo[oo >= a1 & oo <= b1]
	return(result)
}

moveave<-function(v){
	for (i in 3:(length(v)-2)){
		tmp=v[(i-2):(i+2)]
		v[i]=mean(tmp,na.rm=T)
        }
	v[1]=mean(v[1:4],na.rm=T)
	v[2]=mean(v[2:5],na.rm=T)
	v[length(v)-1]=mean(v[(length(v)-4):(length(v)-1)])
	v[length(v)]=mean(v[(length(v)-3):length(v)])
	return (v)
}

supermedian<-function(startend,zscore){#update this function to ignore lowMAPQ regions inside a segment
		if (all(is.na(zscore[startend[1]:startend[2]]))){
			return (median(zscore,na.rm=T))
		}else if (sum(is.na(zscore[startend[1]:startend[2]]))/(startend[2]-startend[1])<.67){
			return (median(zscore[startend[1]:startend[2]],na.rm=T))
		}else{
                        return (median(zscore[startend[1]:startend[2]],na.rm=T)*.5)#It's too strong to use the nonLMAPQ wins if LMAPQ accounts for larger proportion;So cut it to half
		}
		#keptidx=startend[1]:startend[2]
		#keptidx=keptidx[!keptidx %in% rmidx]
		#return (ifelse(length(keptidx)>3,median(zscore[keptidx]),median(zscore)))
	}


compress<-function(vec,step=20){
	init_pos=c(seq(1,length(vec),by=step))
	if (!length(vec) %in% init_pos) init_pos=c(init_pos,length(vec))
	mat=vector()
	for (i in 1:(length(init_pos)-1)){
		mat=rbind(mat,c(init_pos[i],init_pos[i+1]))
	}
	medvec=apply(mat,1,supermedian,zscore=vec)
	return (medvec)
}

cumsumwin<-function(medvec,winnum){
	diff_init=diff(c(rep(median(medvec),winnum),medvec,rep(median(medvec),winnum)),lag=winnum)
	diff_wins=c()
	for (i in 1:(length(diff_init)-winnum)) diff_wins=c(diff_wins,median(diff_init[i:(i+winnum)]))
	return (diff_wins)
}

findpeak<-function(v,step){
	peak=c();valley=c()
	v[is.na(v)]=median(v,na.rm=T)
	trimv=quan(v,.05)
	CIup=median(trimv)+1*sd(trimv)
	CIdown=median(trimv)-1*sd(trimv)
	for (i in 2:(length(v))){
		if (i>=5 & i<=(length(v)-5)){
			if (v[i]>max(v[(i-4):(i-1)]) & max(v[(i+1):(i+4)])<v[i] & v[i]>CIup){
            	peak=c(peak,i)
			}else if (v[i]<min(v[(i-4):(i-1)]) & v[i]<min(v[(i+1):(i+4)]) & v[i]<CIdown){
            	valley=c(valley,i)
            }
        }else if (i<5){
        	if (v[i]>max(v[-i][1:(i+4)]) & v[i]>CIup & !any(1:(i+4) %in% peak)){
            	peak=c(peak,i)
            }else if (v[i]<min(v[-i][1:(i+4)]) & v[i]<CIdown & !any(1:(i+4) %in% valley)){
            	valley=c(valley,i)
            }
        }else if (i>length(v)-5){
        	if (v[i]>max(v[-i][(i-4):(length(v)-1)]) & v[i]>CIup & !any((i-4):length(v) %in% peak)){
            	peak=c(peak,i)
            }else if (v[i]<min(v[-i][(i-4):(length(v)-1)]) & v[i]<CIdown & !any((i-4):length(v) %in% valley)){
                 valley=c(valley,i)
            }
        }
    }

    bp=sort(c(peak,valley))
    if (length(bp)>2){
		rmbp=c()
		for (i in 2:length(bp)){
			if (!all(c(bp[i-1],bp[i+1]) %in% peak) & !all(c(bp[i-1],bp[i+1]) %in% valley)){
				rmbp=c(rmbp,i)
			}
        }
        bp=bp[-rmbp]
    }
    return (sort(unique(c(c(peak-1,valley-1)*step,1,length(v)*step))))#a vector consists of breakpoints such as [1,50,80,190]
}

#PeakSeg<-function(vec,winnum=3,step=20){
    medvec=compress(vec,step)
    medvec=moveave(medvec)
    forward=cumsumwin(medvec,winnum)
    backward=-rev(cumsumwin(rev(medvec),winnum))
    backward=c(rep(NA,5),backward)
    residule=999
    slide=1
    for (i in 1:5){
		res_=sum(abs(forward-backward[i:(i+length(forward)-1)]),na.rm=T)
		if (res_<residule){
			residule=res_
			slide=i
		}
    }
    backward=backward[slide:(slide+length(forward)-1)]
    waves=forward+backward #forward+backward to enhance peak and valley signal
    waves[is.na(waves)]=forward[is.na(waves)]*2
    breakpoints=findpeak(waves,step)
    SegmentsTable=c()
    Segmented=c()
    for (i in 1:(length(breakpoints)-1)){
        from=breakpoints[i]
        to=min(length(vec),breakpoints[i+1])
        SegmentsTable=rbind(SegmentsTable,c(from,to-from+1,median(vec[from:to])))
        Segmented=c(Segmented,rep(median(vec[from:to]),to-from+1))
    }
    seg.data=list(SegmentsTable=SegmentsTable,Segmented=Segmented)
    return (seg.data)
}








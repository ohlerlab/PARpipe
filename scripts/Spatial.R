simulate = function(un, lengths, func){
  ra<-c()
  counter<-c(lengths[1])
  lengths<-append(lengths, 0)
  for (j in 2:length(lengths)){
    if (lengths[j] == lengths[j - 1]){
      counter<-append(counter, lengths[j])
    }
    else{
      if (((counter[1] - 25)/25) > (length(counter) - 1) && length(counter) > 1){
        boo1<-TRUE
        while(boo1){
          counting<-c(func(counter[1]))
          if (counting[1] >= 37.5 || counting[1] <= counter[1] - 37.5){
            boo1<-FALSE
            for (i in 2:length(counter)){
              boo2<-TRUE
              while(boo2){
                boo2<-FALSE
                add<-func(counter[1])
                for (k in 1:length(counting)){
                  if (add < counting[k] - 25){
                    break
                  }else if (add > counting[k] - 25 && add < counting[k] + 25){
                    boo2<-TRUE
                    break
                  }
                }
                if (!boo2){
                  counting<-append(counting, add)
                }
              }
              counting<-sort(counting)
              if (counting[1] <= 37.5 && counting[length(counting)] >= counter[1] - 37.5 && i != length(counter)){
                boo3<-TRUE
                for (k in 2:length(counting)){
                  if (counting[k] > counting[k - 1] + 50){
                    boo3<-FALSE
                    break
                  }
                }
                if (boo3){
                  boo1<-TRUE
                  break
                }
              }
            }
          }
        }
      }else{
        counting<-c()
        for (i in 1:length(counter)){
          counting<-append(counting, func(counter[1]))
        }
      }
      ra<-append(ra, counting)
      counter<-c(lengths[j])
    }
  }
  ra
}

simulateMore = function(un, res, genco, genma){
  li<-c()
  for (j in res){
    if ((genma[j] - 25)/25 > (length(genco[[j]]) - 1)){
      boo1<-TRUE
      while(boo1){
        co<-c(runif(1) * genma[j])
        if (co[1] >= 37.5 || co[1] <= genma[j] - 37.5){
          boo1<-FALSE
          for (i in 2:length(genco[[j]])){
            boo2<-TRUE
            while(boo2){
              boo2<-FALSE
              add<-runif(1) * genma[j]
              for (k in 1:length(co)){
                if (add < co[k] - 25){
                  break
                }else if (add > co[k] - 25 && add < co[k] + 25){
                  boo2<-TRUE
                  break
                }
              }
              if (!boo2){
                co<-append(co, add)
              }
            }
            co<-sort(co)
            if (co[1] <= 37.5 && co[length(co)] >= genma[j] - 37.5 && i != length(genco[[j]])){
              boo3<-TRUE
              for (k in 2:length(co)){
                if (co[k] > co[k - 1] + 50){
                  boo3<-FALSE
                  break
                }
              }
              if (boo3){
                boo1<-TRUE
                break
              }
            }
          }
        }
      }
    }else{
      co<-runif(length(genco[[j]])) * genma[j]
    }
    li<-append(li, co[2] - co[1])
    if (length(co) > 2){
      for (l in 2:(length(co)-1)){
        li<-append(li, min(co[l] - co[l - 1], co[l + 1] - co[l]))
      }
    }
    li<-append(li, co[length(co)] - co[length(co) - 1])
  }
  li
}

simulateMore2 = function(un, res, genco1, genco2, genma, genli1){
  li1<-c()
  for (j in res){
    if ((genma[j,2] - genma[j,1] - 25)/25 > (length(genco1[[j]]) - 1)){
      boo1<-TRUE
      while(boo1){
        co1<-c(runif(1) * (genma[j,2] - genma[j,1]) + genma[j,1])
        if (co1[1] >= 37.5 + genma[j,1] || co1[1] <= genma[j,2] - 37.5){
          boo1<-FALSE
          for (i in 2:length(genco1[[j]])){
            boo2<-TRUE
            while(boo2){
              boo2<-FALSE
              add<-runif(1) * (genma[j,2] - genma[j,1]) + genma[j,1]
              for (k in 1:length(co1)){
                if (add < co1[k] - 25){
                  break
                }else if (add > co1[k] - 25 && add < co1[k] + 25){
                  boo2<-TRUE
                  break
                }
              }
              if (!boo2){
                co1<-append(co1, add)
              }
            }
            co1<-sort(co1)
            if (co1[1] <= 37.5 + genma[j,1] && co1[length(co1)] >= genma[j] - 37.5 && i != length(genco1[[j]])){
              boo3<-TRUE
              for (k in 2:length(co1)){
                if (co1[k] > co1[k - 1] + 50){
                  boo3<-FALSE
                  break
                }
              }
              if (boo3){
                boo1<-TRUE
                break
              }
            }
          }
        }
      }
    }else{
      co1<-runif(length(genco1[[j]])) * (genma[j,2] - genma[j,1]) + genma[j,1]
    }
    for (l in 1:length(co1)){
      li1<-append(li1, -(genco2[[j]] - co1[l])[which.min(sapply(genco2[[j]] - co1[l], abs))])
    }
  }
  li1
}

plotSample = function(points, lengths, dist, title, beg, limit){
  func = function(n) sample((25/2):(n-(25/2)),1)
  if (dist){
    bw<-25
  }else{
    bw<-.005
  }
  r<-range(density(points)$x)
  r[2]<-min(r[2], 10000)
  if (limit > 0){
    bw<-bw/5
    rang<-c(15, limit*(bw/5) + 15)
  }else if (dist == TRUE){
    rang<-c(0, r[2]/(bw/5))
  }else{
    rang<-c(0, 1/(bw/5))
    r<-c(0, 1)
  }
  num<-ceiling((r[2]-r[1])/(bw/5)) + 1
  if (beg == FALSE){
    rang<-rev(rang)
  }
  fr<--3*bw
  if(limit > 0){
    to<-3*bw + 2*limit
  }else if (dist == TRUE){
    to<-floor((3*bw + min(max(lengths), 10000))/(bw/5))*(bw/5)
  }else{
    to<-1 + 3*bw
  }
  num<-ceiling((to-fr)/(bw/5))+1
  dens<-density(points, bw=bw, n=num, from=fr, to=to)
  dom<-range(dens$y)
  cmat<-matrix(ncol=num, nrow=0)
  if (!is.matrix(pmat)){
    pmat<<-matrix(ncol=length(lengths), nrow=0)
    vec<-rep(1, 100)    
    lis<-mclapply(vec, simulate, lengths, func)
    for (i in 1:100){
      pmat<<-rbind(pmat, lis[[i]])
      if (!dist){
        lis[[i]]<-lis[[i]]/lengths
      }
      cmat<-rbind(cmat, density(lis[[i]],bw=bw, n=num, from=fr, to=to)$y)
    }
  }else{
    if (!dist){
      for (i in 1:100){
        pmat[i,]<-pmat[i,]/lengths
        cmat<-rbind(cmat, density(pmat[i,],bw=bw, n=num, from=fr, to=to)$y)
      }
    }else if (!beg){
      for (i in 1:100){
        pmat[i,]<-lengths - pmat[i,]
        cmat<-rbind(cmat, density(pmat[i,],bw=bw, n=num, from=fr, to=to)$y)
      }
    }else{
      for (i in 1:100){
        cmat<-rbind(cmat, density(pmat[i,],bw=bw, n=num, from=fr, to=to)$y)
      }
    }
  }
  clusterplot(cmat, xaxt="n", main=title, xlim=rang, ylim=dom, cex.axis=4/resolution, cex.lab=4/resolution, cex.main=4/resolution, xlab=paste0("N = ", length(points), " Bandwidth = ", bw), ylab="Density", colpal="crazyblue", size=FALSE)
  tick<-axTicks(side=1) + 15
  xp<-tick*(bw/5) + fr
  axis(side=1, labels=xp, at<-tick, cex.axis=4/resolution)
  dens<-density(points,bw=bw, n=num, from=fr, to=to)
  dens$x<-dens$x/(bw/5) + 16
  lines(dens, col="red", lwd=4/resolution)
}

plotMore = function(genli, genco, temp, genma, title, limit){
  if (!is.null(genli)){
    if (length(genli) > 25){
      bw<-25
      r<-range(density(genli)$x)
      r[2]<-min(r[2], 10000)
      if (limit > 0){
        bw<-bw/5
        rang<-c(15, limit*(bw/5) + 15)
      }else{
        rang<-c(0,r[2]/(bw/5))
      }
      num<-ceiling((r[2]-r[1])/(bw/5)) + 1
      fr<--3*bw
      if(limit > 0){
        to<-3*bw + 2*limit
      }else{
        to<-floor((3*bw + max(r))/(bw/5))*(bw/5)
      }
      num<-ceiling((to-fr)/(bw/5))+1
      dens<-density(genli, bw=bw, n=num, from=fr, to=to)
      dom<-range(dens$y)
      cmat<-matrix(ncol=num, nrow=0)
      if (!is.matrix(pmat)){
        res<-which(sapply(genco, length) > 1)
        pmat<<-matrix(ncol=length(genli), nrow=0)
        vec<-rep(1, 100)
        lis<-mclapply(vec, simulateMore, res, genco, genma)
        for (i in 1:100){
          pmat<<-rbind(pmat, lis[[i]])
          cmat<-rbind(cmat, density(lis[[i]], bw=bw, n=num, from=fr, to=to)$y)
        }
      }else{
        for (i in 1:100){
          cmat<-rbind(cmat, density(pmat[i,], bw=bw, n=num, from=fr, to=to)$y)
        }
      }
      clusterplot(cmat, xaxt="n", main=title, xlim=rang, ylim=dom, cex.axis=4/resolution, cex.lab=4/resolution, cex.main=4/resolution, xlab=paste0("N = ", length(genli), " Bandwidth = ", bw), ylab="Density", colpal="crazyblue", size=FALSE)
      tick<-axTicks(side=1) + 15
      xp<-tick*(bw/5) + fr
      axis(side=1, labels=xp, at<-tick, cex.axis=4/resolution)
      dens<-density(genli,bw=bw, n=num, from=fr, to=to)
      dens$x<-dens$x/(bw/5) + 16
      lines(dens, col="red", lwd=4/resolution)
    }
  }
}

plotMore2 = function(genli1, genli2, genco1, genco2, temp1, temp2, genma, title1, title2, limit){
  res1<-which(sapply(genco1, length) > 0)
  res2<-which(sapply(genco2, length) > 0)
  res<-intersect(res1, res2)
  if (length(res) > 0){
    bw<-25
    if (limit > 0){
      bw<-bw/5
    }
    if (length(genli1) > 25){
      r<-range(density(genli1)$x)
      r[2]<-min(r[2], 10000)
      r[1]<-max(r[1], -10000)
      num<-ceiling((r[2]-r[1])/(bw))*5
      if(limit > 0){
        fr<--3*bw - 2*limit
        to<-3*bw + 2*limit
      }else{
        fr<-ceiling((-3*bw + r[1])/(bw/5))*(bw/5)
        to<-floor((3*bw + r[2])/(bw/5))*(bw/5)
      }
      num<-ceiling((to-fr)/(bw/5)) + 1
      if (limit > 0){
        rang<-c(limit*(bw/5) + 15, 3*limit*(bw/5) + 15)
      }else{
        rang<-c(15, num - 15)
      }
      dens<-density(genli1, bw=bw, n=num, from=fr, to=to)
      dom<-range(dens$y)
      cmat<-matrix(ncol=num, nrow=0)
      if (!is.matrix(pmat1)){
        pmat1<<-matrix(ncol=length(genli1), nrow=0)
        vec<-rep(1, 100)
        lis<-mclapply(vec, simulateMore2, res, genco1, genco2, genma, genli1)
        for (i in 1:100){
          pmat1<<-rbind(pmat1, lis[[i]])
          cmat<-rbind(cmat, density(lis[[i]], bw=bw, n=num, from=fr, to=to)$y)
        }
      }else{
        for (i in 1:100){
          cmat<-rbind(cmat, density(pmat1[i,], bw=bw, n=num, from=fr, to=to)$y)
        }
      }
      clusterplot(cmat, xaxt="n", main=title1, xlim=rang, ylim=dom, cex.axis=4/resolution, cex.lab=4/resolution, cex.main=3/resolution, xlab=paste0("N = ", length(genli1), " Bandwidth = ", bw), ylab="Density", colpal="crazyblue", size=FALSE)
      tick<-axTicks(side=1) + 15
      xp<-tick*(bw/5) + fr
      axis(side=1, labels=xp, at<-tick, cex.axis=4/resolution)
      dens<-density(genli1, bw=bw, n=num, from=fr, to=to)
      dens$x<-dens$x/(bw/5) + (-fr/(bw/5)) + 1
      lines(dens, col="red", lwd=4/resolution)
      r<-range(density(genli2)$x)
      r[2]<-min(r[2], 10000)
      r[1]<-max(r[1], -10000)
      num<-ceiling((r[2]-r[1])/(bw))*5
      if (limit > 0){
        fr<--3*bw - 2*limit
        to<-3*bw + 2*limit
      }else{
        fr<-ceiling((-3*bw + r[1])/(bw/5))*(bw/5)
        to<-floor((3*bw + r[2])/(bw/5))*(bw/5)
      }
      num<-ceiling((to-fr)/(bw/5)) + 1
      if (limit > 0){
        rang<-c(limit*(bw/5) + 15, 3*limit*(bw/5) + 15)
      }else{
        rang<-c(15, num - 15)
      }
      dens<-density(genli2, bw=bw, n=num, from=fr, to=to)
      dom<-range(dens$y)
      cmat<-matrix(ncol=num, nrow=0)
      if (!is.matrix(pmat2)){
        pmat2<<-matrix(ncol=length(genli2), nrow=0)
        vec<-rep(1, 100)
        lis<-mclapply(vec, simulateMore2, res, genco2, genco1, genma, genli2)
        for (i in 1:100){
          pmat2<<-rbind(pmat2, lis[[i]])
          cmat<-rbind(cmat, density(lis[[i]], bw=bw, n=num, from=fr, to=to)$y)
        }
      }else{
        for (i in 1:100){
          cmat<-rbind(cmat, density(pmat2[i,], bw=bw, n=num, from=fr, to=to)$y)
        }
      }
      clusterplot(cmat, xaxt="n", main=title2, xlim=rang, ylim=dom, cex.axis=4/resolution, cex.lab=4/resolution, cex.main=3/resolution, xlab=paste0("N = ", length(genli2), " Bandwidth = ", bw), ylab="Density", colpal="crazyred", size=FALSE)
      tick<-axTicks(side=1) + 15
      xp<-tick*(bw/5) + fr
      axis(side=1, labels=xp, at<-tick, cex.axis=4/resolution)
      dens<-density(genli2, bw=bw, n=num, from=fr, to=to)
      dens$x<-dens$x/(bw/5) + (-fr/(bw/5)) + 1
      lines(dens, col="blue", lwd=4/resolution)
    }
  }
}

if (yn == "y"){
  file<-paste0("tmp_", fn, "_", fn2, "_spatial.csv")
}else{
  file<-paste0("tmp_", fn, "_spatial.csv")
}
input<-read.csv(file, header=FALSE, blank.lines.skip=FALSE)
input<-as.matrix(input)
if (yn == "y"){
  temp31<-na.omit(as.numeric(input[3,]))
  temp32<-na.omit(as.numeric(input[4,]))
  temp51<-na.omit(as.numeric(input[5,]))
  temp52<-na.omit(as.numeric(input[6,]))
  tempc1<-na.omit(as.numeric(input[7,]))
  tempc2<-na.omit(as.numeric(input[8,]))
  tempi1<-na.omit(as.numeric(input[9,]))
  tempi2<-na.omit(as.numeric(input[10,]))
  temple1<-na.omit(as.numeric(input[11,]))
  temple2<-na.omit(as.numeric(input[12,]))
  templi1<-na.omit(as.numeric(input[13,]))
  templi2<-na.omit(as.numeric(input[14,]))
  genli31<-na.omit(as.numeric(input[15,]))
  genli32<-na.omit(as.numeric(input[16,]))
  genli51<-na.omit(as.numeric(input[17,]))
  genli52<-na.omit(as.numeric(input[18,]))
  genlic1<-na.omit(as.numeric(input[19,]))
  genlic2<-na.omit(as.numeric(input[20,]))
  genliii1<-na.omit(as.numeric(input[21,]))
  genliii2<-na.omit(as.numeric(input[22,]))
  genlile1<-na.omit(as.numeric(input[23,]))
  genlile2<-na.omit(as.numeric(input[24,]))
  genlili1<-na.omit(as.numeric(input[25,]))
  genlili2<-na.omit(as.numeric(input[26,]))
  genma3<-na.omit(as.numeric(input[27,]))
  genma3<-cbind(genma3, na.omit(as.numeric(input[28,])))
  genma5<-na.omit(as.numeric(input[29,]))
  genma5<-cbind(genma5, na.omit(as.numeric(input[30,])))
  genmac<-na.omit(as.numeric(input[31,]))
  genmac<-cbind(genmac, na.omit(as.numeric(input[32,])))
  genmaii<-na.omit(as.numeric(input[33,]))
  genmaii<-cbind(genmaii, na.omit(as.numeric(input[34,])))
  genmale<-na.omit(as.numeric(input[35,]))
  genmale<-cbind(genmale, na.omit(as.numeric(input[36,])))
  genmali<-na.omit(as.numeric(input[37,]))
  genmali<-cbind(genmali, na.omit(as.numeric(input[38,])))
  genco31 = list();
  tgenco31<-na.omit(input[39,])
  c<-which(tgenco31 == "")
  if (length(c) > 0){
    tgenco31<-tgenco31[-c]
  }
  for (i in tgenco31){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      genco31[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      genco31[[t[1]]]<-numeric(0)
    }
  }
  genco32 = list();
  tgenco32<-na.omit(input[40,])
  c<-which(tgenco32 == "")
  if (length(c) > 0){
    tgenco32<-tgenco32[-c]
  }
  for (i in tgenco32){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      genco32[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      genco32[[t[1]]]<-numeric(0)
    }
  }
  genco51 = list();
  tgenco51<-na.omit(input[41,])
  c<-which(tgenco51 == "")
  if (length(c) > 0){
    tgenco51<-tgenco51[-c]
  }
  for (i in tgenco51){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      genco51[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      genco51[[t[1]]]<-numeric(0)
    }
  }
  genco52 = list();
  tgenco52<-na.omit(input[42,])
  c<-which(tgenco52 == "")
  if (length(c) > 0){
    tgenco52<-tgenco52[-c]
  }
  for (i in tgenco52){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      genco52[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      genco52[[t[1]]]<-numeric(0)
    }
  }
  gencoc1 = list();
  tgencoc1<-na.omit(input[43,])
  c<-which(tgencoc1 == "")
  if (length(c) > 0){
    tgencoc1<-tgencoc1[-c]
  }
  for (i in tgencoc1){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencoc1[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencoc1[[t[1]]]<-numeric(0)
    }
  }
  gencoc2 = list();
  tgencoc2<-na.omit(input[44,])
  c<-which(tgencoc2 == "")
  if (length(c) > 0){
    tgencoc2<-tgencoc2[-c]
  }
  for (i in tgencoc2){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencoc2[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencoc2[[t[1]]]<-numeric(0)
    }
  }
  gencoii1 = list();
  tgencoii1<-na.omit(input[45,])
  c<-which(tgencoii1 == "")
  if (length(c) > 0){
    tgencoii1<-tgencoii1[-c]
  }
  for (i in tgencoii1){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencoii1[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencoii1[[t[1]]]<-numeric(0)
    }
  }
  gencoii2 = list();
  tgencoii2<-na.omit(input[46,])
  c<-which(tgencoii2 == "")
  if (length(c) > 0){
    tgencoii2<-tgencoii2[-c]
  }
  for (i in tgencoii2){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencoii2[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencoii2[[t[1]]]<-numeric(0)
    }
  }
  gencole1 = list();
  tgencole1<-na.omit(input[47,])
  c<-which(tgencole1 == "")
  if (length(c) > 0){
    tgencole1<-tgencole1[-c]
  }
  for (i in tgencole1){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencole1[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencole1[[t[1]]]<-numeric(0)
    }
  }
  gencole2 = list();
  tgencole2<-na.omit(input[48,])
  c<-which(tgencole2 == "")
  if (length(c) > 0){
    tgencole2<-tgencole2[-c]
  }
  for (i in tgencole2){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencole2[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencole2[[t[1]]]<-numeric(0)
    }
  }
  gencoli1 = list();
  tgencoli1<-na.omit(input[49,])
  c<-which(tgencoli1 == "")
  if (length(c) > 0){
    tgencoli1<-tgencoli1[-c]
  }
  for (i in tgencoli1){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencoli1[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencoli1[[t[1]]]<-numeric(0)
    }
  }
  gencoli2 = list();
  tgencoli2<-na.omit(input[50,])
  c<-which(tgencoli2 == "")
  if (length(c) > 0){
    tgencoli2<-tgencoli2[-c]
  }
  for (i in tgencoli2){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencoli2[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencoli2[[t[1]]]<-numeric(0)
    }
  }
  genl31<-na.omit(as.numeric(input[51,]))
  genl32<-na.omit(as.numeric(input[52,]))
  genl51<-na.omit(as.numeric(input[53,]))
  genl52<-na.omit(as.numeric(input[54,]))
  genlc1<-na.omit(as.numeric(input[55,]))
  genlc2<-na.omit(as.numeric(input[56,]))
  genlii1<-na.omit(as.numeric(input[57,]))
  genlii2<-na.omit(as.numeric(input[58,]))
  genlle1<-na.omit(as.numeric(input[59,]))
  genlle2<-na.omit(as.numeric(input[60,]))
  genlli1<-na.omit(as.numeric(input[61,]))
  genlli2<-na.omit(as.numeric(input[62,]))
}else{
  temp3<-na.omit(as.numeric(input[2,]))
  temp5<-na.omit(as.numeric(input[3,]))
  tempc<-na.omit(as.numeric(input[4,]))
  tempi<-na.omit(as.numeric(input[5,]))
  temple<-na.omit(as.numeric(input[6,]))
  templi<-na.omit(as.numeric(input[7,]))
  temp3e<-na.omit(as.numeric(input[8,]))
  temp5e<-na.omit(as.numeric(input[9,]))
  tempce<-na.omit(as.numeric(input[10,]))
  tempid<-na.omit(as.numeric(input[11,]))
  templed<-na.omit(as.numeric(input[12,]))
  templid<-na.omit(as.numeric(input[13,]))
  genli3<-na.omit(as.numeric(input[14,]))
  genli5<-na.omit(as.numeric(input[15,]))
  genlic<-na.omit(as.numeric(input[16,]))
  genliii<-na.omit(as.numeric(input[17,]))
  genlile<-na.omit(as.numeric(input[18,]))
  genlili<-na.omit(as.numeric(input[19,]))
  genma3<-na.omit(as.numeric(input[20,]))
  genma5<-na.omit(as.numeric(input[21,]))
  genmac<-na.omit(as.numeric(input[22,]))
  genmaii<-na.omit(as.numeric(input[23,]))
  genmale<-na.omit(as.numeric(input[24,]))
  genmali<-na.omit(as.numeric(input[25,]))
  genco3 = list();
  tgenco3<-na.omit(input[26,])
  c<-which(tgenco3 == "")
  if (length(c) > 0){
    tgenco3<-tgenco3[-c]
  }
  for (i in tgenco3){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      genco3[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      genco3[[t[1]]]<-c()
    }
  }
  genco5 = list();
  tgenco5<-na.omit(input[27,])
  c<-which(tgenco5 == "")
  if (length(c) > 0){
    tgenco5<-tgenco5[-c]
  }
  for (i in tgenco5){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      genco5[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      genco5[[t[1]]]<-c()
    }
  }
  gencoc = list();
  tgencoc<-na.omit(input[28,])
  c<-which(tgencoc == "")
  if (length(c) > 0){
    tgencoc<-tgencoc[-c]
  }
  for (i in tgencoc){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencoc[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencoc[[t[1]]]<-c()
    }
  }
  gencoii = list();
  tgencoii<-na.omit(input[29,])
  c<-which(tgencoii == "")
  if(length(c) > 0){
    tgencoii<-tgencoii[-c]
  }
  for (i in tgencoii){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencoii[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencoii[[t[1]]]<-c()
    }
  }
  gencole = list();
  tgencole<-na.omit(input[30,])
  c<-which(tgencole == "")
  if (length(c) > 0){
    tgencole<-tgencole[-c]
  }
  for (i in tgencole){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencole[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencole[[t[1]]]<-c()
    }
  }
  gencoli = list();
  tgencoli<-na.omit(input[31,])
  c<-which(tgencoli == "")
  if (length(c) > 0){
    tgencoli<-tgencoli[-c]
  }
  for (i in tgencoli){
    t<-unlist(strsplit(i, ":"))
    if (length(t) > 1){
      gencoli[[t[1]]]<-as.numeric(t[2:length(t)])
    }else{
      gencoli[[t[1]]]<-c()
    }
  }
}
if (yn != "y"){
  beg<-TRUE
  if (length(temp3) > 25){
    pdf(paste0(fn, "_Spatial_3'utr.pdf"), width=19.25/resolution, height=12.375/resolution)
    par(mar=c(10/resolution, 10/resolution, 10/resolution, 10/resolution), mgp=c(7/resolution, 2.5/resolution, 0))
    pmat<-c()
    plotSample(temp3, temp3e, TRUE, "Distance Downstream of Stop Codon", beg, -1)
    plotSample(temp3, temp3e, TRUE, "Distance Downstream of Stop Codon", beg, 200)
    temp3<-temp3e-temp3
    plotSample(temp3, temp3e, TRUE, "Distance Upstream of 3'UTR End", !beg, -1)
    plotSample(temp3, temp3e, TRUE, "Distance Upstream of 3'UTR End", !beg, 200)
    temp3<-(temp3e-temp3)/temp3e
    plotSample(temp3, temp3e, FALSE, "Normalized Distribution of Sites: 3'UTR", beg, -1)
    pmat<-c()
    plotMore(genli3, genco3, temp3, genma3, "Distance Between Sites in the 3'UTR", -1)
    plotMore(genli3, genco3, temp3, genma3, "Distance Between Sites in the 3'UTR", 100)
    dev.off()
  }
  if (length(temp5) > 25){
    pdf(paste0(fn, "_Spatial_5'utr.pdf"), width=19.25/resolution, height=12.375/resolution)
    par(mar=c(10/resolution, 10/resolution, 10/resolution, 10/resolution), mgp=c(7/resolution, 2.5/resolution, 0))
    pmat<-c()
    plotSample(temp5, temp5e, TRUE, "Distance Downstream of 5'UTR Beginning", beg, -1)
    plotSample(temp5, temp5e, TRUE, "Distance Downstream of 5'UTR Beginning", beg, 200)
    temp5<-temp5e-temp5
    plotSample(temp5, temp5e, TRUE, "Distance Upstream of Start Codon", !beg, -1)
    plotSample(temp5, temp5e, TRUE, "Distance Upstream of Start Codon", !beg, 200)
    temp5<-(temp5e-temp5)/temp5e
    plotSample(temp5, temp5e, FALSE, "Normalized Distribution of Sites: 5'UTR", beg, -1)
    pmat<-c()
    plotMore(genli5, genco5, temp5, genma5, "Distance Between Sites in the 5'UTR", -1)
    plotMore(genli5, genco5, temp5, genma5, "Distance Between Sites in the 5'UTR", 100)
    dev.off()
  }
  if (length(tempc) > 25){
    pdf(paste0(fn, "_Spatial_Coding.pdf"), width=19.25/resolution, height=12.375/resolution)
    par(mar=c(10/resolution, 10/resolution, 10/resolution, 10/resolution), mgp=c(7/resolution, 2.5/resolution, 0))
    pmat<-c()
    plotSample(tempc, tempce, TRUE, "Distance Downstream of CDS Exon Beginning", beg, -1)
    plotSample(tempc, tempce, TRUE, "Distance Downstream of CDS Exon Beginning", beg, 200)
    tempc<-tempce-tempc
    plotSample(tempc, tempce, TRUE, "Distance Upstream of CDS Exon End", !beg, -1)
    plotSample(tempc, tempce, TRUE, "Distance Upstream of CDS Exon End", !beg, 200)
    tempc<-(tempce-tempc)/tempce
    plotSample(tempc, tempce, FALSE, "Normalized Distribution of Sites: CDS Exon", beg, -1)
    pmat<-c()
    plotMore(genlic, gencoc, tempc, genmac, "Distance Between Sites in the CDS", -1)
    plotMore(genlic, gencoc, tempc, genmac, "Distance Between Sites in the CDS", 100)
    dev.off()
  }
  if (length(tempi) > 25){
    pdf(paste0(fn, "_Spatial_Intron.pdf"), width=19.25/resolution, height=12.375/resolution)
    par(mar=c(10/resolution, 10/resolution, 10/resolution, 10/resolution), mgp=c(7/resolution, 2.5/resolution, 0))
    pmat<-c()
    plotSample(tempi, tempid, TRUE, "Distance Downstream of 5' Splice Site", beg, -1)
    plotSample(tempi, tempid, TRUE, "Distance Downstream of 5' Splice Site", beg, 500)
    tempi<-tempid-tempi
    plotSample(tempi, tempid, TRUE, "Distance Upstream of 3' Splice Site", !beg, -1)
    plotSample(tempi, tempid, TRUE, "Distance Upstream of 3' Splice Site", !beg, 500)
    tempi<-(tempid-tempi)/tempid
    plotSample(tempi, tempid, FALSE, "Normalized Distribution of Sites: Intron", beg, -1)
    pmat<-c()
    plotMore(genliii, gencoii, tempi, genmaii, "Distance Between Sites in the Intron", -1)
    plotMore(genliii, gencoii, tempi, genmaii, "Distance Between Sites in the Intron", 100)
    dev.off()
  }
  if (length(temple) > 25){
    pdf(paste0(fn, "_Spatial_lincRNA_Exon.pdf"), width=19.25/resolution, height=12.375/resolution)
    par(mar=c(10/resolution, 10/resolution, 10/resolution, 10/resolution), mgp=c(7/resolution, 2.5/resolution, 0))
    pmat<-c()
    plotSample(temple, templed, TRUE, "Distance Downstream of lincRNA Exon Beginning", beg, -1)
    plotSample(temple, templed, TRUE, "Distance Downstream of lincRNA Exon Beginning", beg, 200)
    temple<-templed-temple
    plotSample(temple, templed, TRUE, "Distance Upstream of lincRNA Exon End", !beg, -1)
    plotSample(temple, templed, TRUE, "Distance Upstream of lincRNA Exon End", !beg, 200)
    temple<-(templed-temple)/templed
    plotSample(temple, templed, FALSE, "Normalized Distribution of Sites: lincRNA Exon", beg, -1)
    pmat<-c()
    plotMore(genlile, gencole, temple, genmale, "Distance Between Sites in lincRNA Exon", -1)
    plotMore(genlile, gencole, temple, genmale, "Distance Between Sites in lincRNA Exon", 100)
    dev.off()
  }
  if (length(templi) > 25){
    pdf(paste0(fn, "_Spatial_lincRNA_Intron.pdf"), width=19.25/resolution, height=12.375/resolution)
    par(mar=c(10/resolution, 10/resolution, 10/resolution, 10/resolution), mgp=c(7/resolution, 2.5/resolution, 0))
    pmat<-c()
    plotSample(templi, templid, TRUE, "Distance Downstream of lincRNA Intron 5' Splice Site", beg, -1)
    plotSample(templi, templid, TRUE, "Distance Downstream of lincRNA Intron 5' Splice Site", beg, 500)
    templi<-templid-templi
    plotSample(templi, templid, TRUE, "Distance Upstream of lincRNA Intron 3' Splice Site", !beg, -1)
    plotSample(templi, templid, TRUE, "Distance Upstream of lincRNA Intron 3' Splice Site", !beg, 500)
    templi<-(templid-templi)/templid
    plotSample(templi, templid, FALSE, "Normalized Distribution of Sites: lincRNA Intron", beg, -1)
    pmat<-c()
    plotMore(genlili, gencoli, templi, genmali, "Distance Between Sites in lincRNA Intron", -1)
    plotMore(genlili, gencoli, templi, genmali, "Distance Between Sites in lincRNA Intron", 100)
    dev.off()
  }
}else{
  pdf(paste0(fn, "_", fn2, "_Spatial.pdf"), width=19.25/resolution, height=12.375/resolution)
  par(mar=c(10/resolution, 10/resolution, 10/resolution, 10/resolution), mgp=c(7/resolution, 2.5/resolution, 0))
  pmat1<-c()
  pmat2<-c()
  plotMore2(genli31, genli32, genco31, genco32, temp31, temp32, genma3, paste0("Distance ", fn, " vs. ", fn2, ": 3'UTR"), paste0("Distance ", fn2, " vs. ", fn, ": 3'UTR"), -1)
  plotMore2(genli31, genli32, genco31, genco32, temp31, temp32, genma3, paste0("Distance ", fn, " vs. ", fn2, ": 3'UTR"), paste0("Distance ", fn2, " vs. ", fn, ": 3'UTR"), 100)
  pmat1<-c()
  pmat2<-c()
  plotMore2(genli51, genli52, genco51, genco52, temp51, temp52, genma5, paste0("Distance ", fn, " vs. ", fn2, ": 5'UTR"), paste0("Distance ", fn2, " vs. ", fn, ": 5'UTR"), -1)
  plotMore2(genli51, genli52, genco51, genco52, temp51, temp52, genma5, paste0("Distance ", fn, " vs. ", fn2, ": 5'UTR"), paste0("Distance ", fn2, " vs. ", fn, ": 5'UTR"), 100)
  pmat1<-c()
  pmat2<-c()
  plotMore2(genlic1, genlic2, gencoc1, gencoc2, tempc1, tempc2, genmac, paste0("Distance ", fn, " vs. ", fn2, ": Coding"), paste0("Distance ", fn2, " vs. ", fn, ": Coding"), -1)
  plotMore2(genlic1, genlic2, gencoc1, gencoc2, tempc1, tempc2, genmac, paste0("Distance ", fn, " vs. ", fn2, ": Coding"), paste0("Distance ", fn2, " vs. ", fn, ": Coding"), 100)
  pmat1<-c()
  pmat2<-c()
  plotMore2(genliii1, genliii2, gencoii1, gencoii2, tempii1, tempii2, genmaii, paste0("Distance ", fn, " vs. ", fn2, ": Intron"), paste0("Distance ", fn2, " vs. ", fn, ": Intron"), -1)
  plotMore2(genliii1, genliii2, gencoii1, gencoii2, tempii1, tempii2, genmaii, paste0("Distance ", fn, " vs. ", fn2, ": Intron"), paste0("Distance ", fn2, " vs. ", fn, ": Intron"), 100)
  pmat1<-c()
  pmat2<-c()
  plotMore2(genlile1, genlile2, gencole1, gencole2, temple1, temple2, genmale, paste0("Distance ", fn, " vs. ", fn2, ": lincRNA Exon"), paste0("Distance ", fn2, " vs. ", fn, ": lincRNA Exon"), -1)
  plotMore2(genlile1, genlile2, gencole1, gencole2, temple1, temple2, genmale, paste0("Distance ", fn, " vs. ", fn2, ": lincRNA Exon"), paste0("Distance ", fn2, " vs. ", fn, ": lincRNA Exon"), 100)
  pmat1<-c()
  pmat2<-c()
  plotMore2(genlili1, genlili2, gencoli1, gencoli2, templi1, templi2, genmali, paste0("Distance ", fn, " vs. ", fn2, ": lincRNA Intron"), paste0("Distance ", fn2, " vs. ", fn, ": lincRNA Intron"), -1)
  plotMore2(genlili1, genlili2, gencoli1, gencoli2, templi1, templi2, genmali, paste0("Distance ", fn, " vs. ", fn2, ": lincRNA Intron"), paste0("Distance ", fn2, " vs. ", fn, ": lincRNA Intron"), 100)
  dev.off()
}

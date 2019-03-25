# This file was produced by S. Guindon.

# Pre-processing of output files

```javascript
cd /home/guindon/latex/student/adrien/revision/res/merged_2018_12_26/AA-900
cat Arbre*/Output_merged_nodes.txt > all
perl -p -i -e 's/\t/\t\t/g' all
perl -p -i -e 's/ /\t\t/g' all
perl -p -i -e 's/\|//g' all
perl -p -i -e 's/Nodes\t\tID/# ID /g' all
perl -p -i -e 's/\//\t\t/g' all
perl -p -i -e 's/\-/\t\t/g' all
Use emacs to remove last four columns
perl -p -i -e 's/\t\t/ /g' all

Remove # from first row.

cd /home/guindon/latex/student/adrien/revision/res/merged_2018_12_26/Codon-900-withProba
cat Arbre*/Output_merged_prob.txt > all
perl -p -i -e 's/\t/\t\t/g' all
perl -p -i -e 's/ /\t\t/g' all
perl -p -i -e 's/\|//g' all
perl -p -i -e 's/Nodes\t\tID/# ID /g' all
perl -p -i -e 's/\//\t\t/g' all
perl -p -i -e 's/\-/\t\t/g' all
Use emacs to remove last four columns
perl -p -i -e 's/\t\t/ /g' all
```

Remove # from first row.

```{r}
library(pbapply);

source("/home/guindon/latex/student/adrien/revision/res/ancestral_analyse.R");
filename="/home/guindon/latex/student/adrien/revision/res/merged_2018_12_26/Codon-900-withProba/all";
d=read.table(filename,header=T,comment="#");
p = as.matrix(cbind(d$pA,d$pC,d$pG,d$pT));
s = pbapply(p,1,mpee,mesh.size=50);
d$mpee=s;
s = pbapply(p,1,which.map);
d$map=s;
s = pbapply(p,1,brier);
d$brier=s;
s = pbapply(p,1,map.plus,1/4);
d$mapplus=s;
s = pbapply(p,1,cum.prob,0.9);
d$cumprob=s;
s = pbapply(p,1,prob.thresh,1/4);
d$probthresh=s;
write.table(d,file=paste(filename,"_checked",sep=""));

source("/home/guindon/latex/student/adrien/revision/res/ancestral_analyse.R");
filename="/home/guindon/latex/student/adrien/revision/res/merged_2018_12_26/AA-900/all";
d=read.table(filename,header=T,comment="#");
p = as.matrix(cbind(d$pA,d$pR,d$pN,d$pD,d$pC,d$pQ,d$pE,d$pG,d$pH,d$pI,d$pL,d$pK,d$pM,d$pF,d$pP,d$pS,d$pT,d$pW,d$pY,d$pV));
s = pbapply(p,1,mpee,mesh.size=50);
d$mpee=s;
s = pbapply(p,1,which.map);
d$map=s;
s = pbapply(p,1,brier);
d$brier=s;
s = pbapply(p,1,map.plus,1/20);
d$mapplus=s;
s = pbapply(p,1,cum.prob,0.9);
d$cumprob=s;
s = pbapply(p,1,prob.thresh,1/20);
d$probthresh=s;
write.table(d,file=paste(filename,"_checked",sep=""));










# Plots

source("/home/guindon/latex/student/adrien/revision/res/ancestral_analyse.R");
## d = read.table("~/latex/student/adrien/revision/res/merged_2018_12_26/Codon-900-withProba/all_checked",header=T);
## d = d[apply(cbind(d$pA,d$pC,d$pG,d$pT),1,max) < 0.95,];
d = read.table("~/latex/student/adrien/revision/res/merged_2018_12_26/AA-900/all_checked",header=T);
d = d[apply(cbind(d$pA,d$pR,d$pN,d$pD,d$pC,d$pQ,d$pE,d$pG,d$pH,d$pI,d$pL,d$pK,d$pM,d$pF,d$pP,d$pS,d$pT,d$pW,d$pY,d$pV),1,max) < 0.95,];

mpee.accu = accuracy(as.vector(d$True),as.vector(d$mpee));
mpee.prec = precision(as.vector(d$mpee));

brier.accu = accuracy(as.vector(d$True),as.vector(d$brier));
brier.prec = precision(as.vector(d$brier));

mapplus.accu = accuracy(as.vector(d$True),as.vector(d$mapplus));
mapplus.prec = precision(as.vector(d$mapplus));

cumprob.accu = accuracy(as.vector(d$True),as.vector(d$cumprob));
cumprob.prec = precision(as.vector(d$cumprob));

probthresh.accu = accuracy(as.vector(d$True),as.vector(d$probthresh));
probthresh.prec = precision(as.vector(d$probthresh));


map.accu = accuracy(as.vector(d$True),as.vector(d$map));
map.prec = precision(as.vector(d$map));

par(mfrow=c(2,4));

col.mpee=rgb(0.6,0.4,0.1,0.5);
col.brier=rgb(0.1,0.6,0.3,0.5);
col.mapplus=rgb(0.4,0.2,0.6,0.5);
col.cumprob=rgb(0.4,0.5,0.2,0.5);
col.probthresh=rgb(0.1,0.3,0.7,0.4);
col.map=rgb(0.2,0.2,0.2,0.5);



amb = 1;
a11 = sum(mpee.accu[mpee.prec == amb] == 1);
b11 = sum(brier.accu[brier.prec == amb] == 1);
c11 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d11 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e11 = sum(probthresh.accu[probthresh.prec == amb] == 1);
f11 = sum(map.accu[map.prec == amb] == 1);

a01 = sum(mpee.accu[mpee.prec == amb] == 0);
b01 = sum(brier.accu[brier.prec == amb] == 0);
c01 = sum(mapplus.accu[mapplus.prec == amb] == 0);
d01 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e01 = sum(probthresh.accu[probthresh.prec == amb] == 0);
f01 = sum(map.accu[map.prec == amb] == 0);


amb = 2;
a12 = sum(mpee.accu[mpee.prec == amb] == 1);
b12 = sum(brier.accu[brier.prec == amb] == 1);
c12 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d12 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e12 = sum(probthresh.accu[probthresh.prec == amb] == 1);
f12 = sum(map.accu[map.prec == amb] == 1);

a02 = sum(mpee.accu[mpee.prec == amb] == 0);
b02 = sum(brier.accu[brier.prec == amb] == 0);
c02 = sum(mapplus.accu[mapplus.prec == amb] == 0);
d02 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e02 = sum(probthresh.accu[probthresh.prec == amb] == 0);
f02 = sum(map.accu[map.prec == amb] == 0);



amb = 3;
a13 = sum(mpee.accu[mpee.prec == amb] == 1);
b13 = sum(brier.accu[brier.prec == amb] == 1);
c13 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d13 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e13 = sum(probthresh.accu[probthresh.prec == amb] == 1);
f13 = sum(map.accu[map.prec == amb] == 1);

a03 = sum(mpee.accu[mpee.prec == amb] == 0);
b03 = sum(brier.accu[brier.prec == amb] == 0);
c03 = sum(mapplus.accu[mapplus.prec == amb] == 0);
d03 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e03 = sum(probthresh.accu[probthresh.prec == amb] == 0);
f03 = sum(map.accu[map.prec == amb] == 0);



amb = 4;
a14 = sum(mpee.accu[mpee.prec >= amb] == 1);
b14 = sum(brier.accu[brier.prec >= amb] == 1);
c14 = sum(mapplus.accu[mapplus.prec >= amb] == 1);
d14 = sum(cumprob.accu[cumprob.prec >= amb] == 1);
e14 = sum(probthresh.accu[probthresh.prec >= amb] == 1);
f14 = sum(map.accu[map.prec == amb] >= 1);

a04 = sum(mpee.accu[mpee.prec >= amb] == 0);
b04 = sum(brier.accu[brier.prec >= amb] == 0);
c04 = sum(mapplus.accu[mapplus.prec >= amb] == 0);
d04 = sum(cumprob.accu[cumprob.prec >= amb] == 0);
e04 = sum(probthresh.accu[probthresh.prec >= amb] == 0);
f04 = sum(map.accu[map.prec >= amb] == 0);



cat(a01," ",a11," ",a01/(a11+a01)," ",(a11+a01)/(a11+a01+a12+a02+a13+a03+a14+a04),"\n");
cat(b01," ",b11," ",b01/(b11+b01)," ",(b11+b01)/(b11+b01+b12+b02+b13+b03+b14+b04),"\n");
cat(c01," ",c11," ",c01/(c11+c01)," ",(c11+c01)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d01," ",d11," ",d01/(d11+d01)," ",(d11+d01)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e01," ",e11," ",e01/(e11+e01)," ",(e11+e01)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");
cat(f01," ",f11," ",f01/(f11+f01)," ",(f11+f01)/(f11+f01+f12+f02+f13+f03+f14+f04),"\n");


cat(a02," ",a12," ",a02/(a12+a02)," ",(a12+a02)/(a11+a01+a12+a02+a13+a03+a14+a04),"\n");
cat(b02," ",b12," ",b02/(b12+b02)," ",(b12+b02)/(b11+b01+b12+b02+b13+b03+b14+b04),"\n");
cat(c02," ",c12," ",c02/(c12+c02)," ",(c12+c02)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d02," ",d12," ",d02/(d12+d02)," ",(d12+d02)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e02," ",e12," ",e02/(e12+e02)," ",(e12+e02)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");
cat(f02," ",f12," ",f02/(f12+f02)," ",(f12+f02)/(f11+f01+f12+f02+f13+f03+f14+f04),"\n");

cat(a03," ",a13," ",a03/(a13+a03)," ",(a13+a03)/(a11+a01+a12+a02+a13+a03+a14+a04),"\n");
cat(b03," ",b13," ",b03/(b13+b03)," ",(b13+b03)/(b11+b01+b12+b02+b13+b03+b14+b04),"\n");
cat(c03," ",c13," ",c03/(c13+c03)," ",(c13+c03)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d03," ",d13," ",d03/(d13+d03)," ",(d13+d03)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e03," ",e13," ",e03/(e13+e03)," ",(e13+e03)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");
cat(f03," ",f13," ",f03/(f13+f03)," ",(f13+f03)/(f11+f01+f12+f02+f13+f03+f14+f04),"\n");

cat(a04," ",a14," ",a04/(a14+a04)," ",(a14+a04)/(a11+a01+a12+a02+a13+a03+a14+a04),"\n");
cat(b04," ",b14," ",b04/(b14+b04)," ",(b14+b04)/(b11+b01+b12+b02+b13+b03+b14+b04),"\n");
cat(c04," ",c14," ",c04/(c14+c04)," ",(c14+c04)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d04," ",d14," ",d04/(d14+d04)," ",(d14+d04)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e04," ",e14," ",e04/(e14+e04)," ",(e14+e04)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");
cat(f04," ",f14," ",f04/(f14+f04)," ",(f14+f04)/(f11+f01+f12+f02+f13+f03+f14+f04),"\n");


barplot(-c(a11,b11,c11,d11,e11,f11),horiz=T,col=c(col.mpee,col.brier,col.mapplus,col.cumprob,col.probthresh,col.map),xlim=c(-max(a11,b11,c11,d11,e11,f11,a01,b01,c01,d01,e01,f01),0),space=0);
barplot(c(a01,b01,c01,d01,e01,f01),horiz=T,col=c(col.mpee,col.brier,col.mapplus,col.cumprob,col.probthresh,col.map),xlim=c(0,max(a11,b11,c11,d11,e11,f11,a01,b01,c01,d01,e01,f01)),space=0);

barplot(-c(a12,b12,c12,d12,e12,f12),horiz=T,col=c(col.mpee,col.brier,col.mapplus,col.cumprob,col.probthresh,col.map),xlim=c(-max(a12,b12,c12,d12,e12,f12,a02,b02,c02,d02,e02,f02),0),space=0);
barplot(c(a02,b02,c02,d02,e02,f02),horiz=T,col=c(col.mpee,col.brier,col.mapplus,col.cumprob,col.probthresh,col.map),xlim=c(0,max(a12,b12,c12,d12,e12,f12,a02,b02,c02,d02,e02,f02)),space=0);

barplot(-c(a13,b13,c13,d13,e13,f13),horiz=T,col=c(col.mpee,col.brier,col.mapplus,col.cumprob,col.probthresh,col.map),xlim=c(-max(a13,b13,c13,d13,e13,f13,a03,b03,c03,d03,e03,f03),0),space=0);
barplot(c(a03,b03,c03,d03,e03,f03),horiz=T,col=c(col.mpee,col.brier,col.mapplus,col.cumprob,col.probthresh,col.map),xlim=c(0,max(a13,b13,c13,d13,e13,f13,a03,b03,c03,d03,e03,f03)),space=0);

barplot(-c(a14,b14,c14,d14,e14,f14),horiz=T,col=c(col.mpee,col.brier,col.mapplus,col.cumprob,col.probthresh,col.map),xlim=c(-max(a14,b14,c14,d14,e14,f14,a04,b04,c04,d04,e04,f04),0),space=0);
barplot(c(a04,b04,c04,d04,e04,f04),horiz=T,col=c(col.mpee,col.brier,col.mapplus,col.cumprob,col.probthresh,col.map),xlim=c(0,max(a14,b14,c14,d14,e14,f14,a04,b04,c04,d04,e04,f04)),space=0);





# Extra results


source("/home/guindon/latex/student/adrien/revision/res/ancestral_analyse.R");
d = read.table("~/latex/student/adrien/revision/res/merged_2018_12_26/Codon-900-withProba/all",header=T);
d = d[apply(cbind(d$pA,d$pC,d$pG,d$pT),1,max) < 0.95,];
p = as.matrix(cbind(d$pA,d$pC,d$pG,d$pT));

diff = pbapply(p,1,map.plus,0.05);
cumprob=pbapply(p,1,cum.prob,0.5);
thresh=pbapply(p,1,prob.thresh,0.45);

diff.accu=accuracy(as.vector(d$True),as.vector(diff));
diff.prec=precision(as.vector(diff));


cumprob.accu=accuracy(as.vector(d$True),as.vector(cumprob));
cumprob.prec=precision(as.vector(cumprob));

thresh.accu=accuracy(as.vector(d$True),as.vector(thresh));
thresh.prec=precision(as.vector(thresh));


amb = 1;
c11 = sum(diff.accu[diff.prec == amb] == 1);
d11 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e11 = sum(thresh.accu[thresh.prec == amb] == 1);

c01 = sum(diff.accu[diff.prec == amb] == 0);
d01 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e01 = sum(thresh.accu[thresh.prec == amb] == 0);


amb = 2;
c12 = sum(diff.accu[diff.prec == amb] == 1);
d12 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e12 = sum(thresh.accu[thresh.prec == amb] == 1);

c02 = sum(diff.accu[diff.prec == amb] == 0);
d02 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e02 = sum(thresh.accu[thresh.prec == amb] == 0);



amb = 3;
c13 = sum(diff.accu[diff.prec == amb] == 1);
d13 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e13 = sum(thresh.accu[thresh.prec == amb] == 1);

c03 = sum(diff.accu[diff.prec == amb] == 0);
d03 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e03 = sum(thresh.accu[thresh.prec == amb] == 0);



amb = 4;
c14 = sum(diff.accu[diff.prec >= amb] == 1);
d14 = sum(cumprob.accu[cumprob.prec >= amb] == 1);
e14 = sum(thresh.accu[thresh.prec >= amb] == 1);

c04 = sum(diff.accu[diff.prec >= amb] == 0);
d04 = sum(cumprob.accu[cumprob.prec >= amb] == 0);
e04 = sum(thresh.accu[thresh.prec >= amb] == 0);



cat(c01," ",c11," ",c01/(c11+c01)," ",(c11+c01)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d01," ",d11," ",d01/(d11+d01)," ",(d11+d01)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e01," ",e11," ",e01/(e11+e01)," ",(e11+e01)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");


cat(c02," ",c12," ",c02/(c12+c02)," ",(c12+c02)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d02," ",d12," ",d02/(d12+d02)," ",(d12+d02)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e02," ",e12," ",e02/(e12+e02)," ",(e12+e02)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c03," ",c13," ",c03/(c13+c03)," ",(c13+c03)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d03," ",d13," ",d03/(d13+d03)," ",(d13+d03)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e03," ",e13," ",e03/(e13+e03)," ",(e13+e03)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c04," ",c14," ",c04/(c14+c04)," ",(c14+c04)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d04," ",d14," ",d04/(d14+d04)," ",(d14+d04)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e04," ",e14," ",e04/(e14+e04)," ",(e14+e04)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");



# RAxML

## DNA
d = read.table("~/latex/student/adrien/revision/res/merged_2018_12_26/Codon-900-withProba/all",header=T);
d = d[apply(cbind(d$pA,d$pC,d$pG,d$pT),1,max) < 0.95,];
raxml.accu=accuracy(as.vector(d$True),as.vector(d$RAxML));
raxml.prec=precision(as.vector(d$RAxML));
sum(raxml.accu[raxml.prec == 1] == 0)/(sum(raxml.accu[raxml.prec == 1] == 0) + sum(raxml.accu[raxml.prec == 1] == 1));

## AA
d = read.table("~/latex/student/adrien/revision/res/merged_2018_12_26/AA-900/all",header=T);
d = d[apply(cbind(d$pA,d$pR,d$pN,d$pD,d$pC,d$pQ,d$pE,d$pG,d$pH,d$pI,d$pL,d$pK,d$pM,d$pF,d$pP,d$pS,d$pT,d$pW,d$pY,d$pV),1,max) < 0.95,];
raxml.accu=accuracy(as.vector(d$True),as.vector(d$RAxML));
raxml.prec=precision(as.vector(d$RAxML));
sum(raxml.accu[raxml.prec == 1] == 0)/(sum(raxml.accu[raxml.prec == 1] == 0) + sum(raxml.accu[raxml.prec == 1] == 1));




# Brier scores

source("/home/guindon/latex/student/adrien/revision2/res/ancestral_analyse.R");

d = read.table("~/latex/student/adrien/revision2/res/merged_2018_12_26/Codon-900-withProba/all_checked",header=T);
d = d[apply(cbind(d$pA,d$pC,d$pG,d$pT),1,max) < 0.95,];
p = matrix(1,length(d$pA),4);
mpee.brier.nt = brier.score.all(as.vector(d$True),as.vector(d$mpee),p,4);
brier.brier.nt = brier.score.all(as.vector(d$True),as.vector(d$brier),p,4);
mapplus.brier.nt = brier.score.all(as.vector(d$True),as.vector(d$mapplus),p,4);
cumprob.brier.nt = brier.score.all(as.vector(d$True),as.vector(d$cumprob),p,4);
probthresh.brier.nt = brier.score.all(as.vector(d$True),as.vector(d$probthresh),p,4);
map.brier.nt = brier.score.all(as.vector(d$True),as.vector(d$map),p,4);
raxml.brier.nt = brier.score.all(as.vector(d$True),as.vector(d$RAxML),p,4);
rand.brier.nt = brier.score.all(as.vector(d$True),sample(c("A","C","G","T"),length(d$True),replace=TRUE),p,4);
p = cbind(d$pA,d$pC,d$pG,d$pT);
post.brier.nt = brier.score.all(as.vector(d$True),as.vector(rep("ACGT",length(d$True))),p,4);

mean(mpee.brier.nt);
mean(brier.brier.nt);
mean(cumprob.brier.nt);
mean(probthresh.brier.nt);
mean(mapplus.brier.nt);
mean(map.brier.nt);
mean(raxml.brier.nt);
mean(post.brier.nt);
mean(rand.brier.nt);


d = read.table("~/latex/student/adrien/revision2/res/merged_2018_12_26/AA-900/all_checked",header=T);
d = d[apply(cbind(d$pA,d$pR,d$pN,d$pD,d$pC,d$pQ,d$pE,d$pG,d$pH,d$pI,d$pL,d$pK,d$pM,d$pF,d$pP,d$pS,d$pT,d$pW,d$pY,d$pV),1,max) < 0.95,];
p = matrix(1,length(d$pA),20);
mpee.brier.aa = brier.score.all(as.vector(d$True),as.vector(d$mpee),p,20);
brier.brier.aa = brier.score.all(as.vector(d$True),as.vector(d$brier),p,20);
mapplus.brier.aa = brier.score.all(as.vector(d$True),as.vector(d$mapplus),p,20);
cumprob.brier.aa = brier.score.all(as.vector(d$True),as.vector(d$cumprob),p,20);
probthresh.brier.aa = brier.score.all(as.vector(d$True),as.vector(d$probthresh),p,20);
map.brier.aa = brier.score.all(as.vector(d$True),as.vector(d$map),p,20);
raxml.brier.aa = brier.score.all(as.vector(d$True),as.vector(d$RAxML),p,20);
rand.brier.aa = brier.score.all(as.vector(d$True),sample(c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"),length(d$True),replace=TRUE),p,20);
p = cbind(d$pA,d$pR,d$pN,d$pD,d$pC,d$pQ,d$pE,d$pG,d$pH,d$pI,d$pL,d$pK,d$pM,d$pF,d$pP,d$pS,d$pT,d$pW,d$pY,d$pV);
post.brier.aa = brier.score.all(as.vector(d$True),as.vector(rep("ARNDCQEGHILKMFPSTWYV",length(d$True))),p,20);


mean(mpee.brier.aa);
mean(brier.brier.aa);
mean(cumprob.brier.aa);
mean(probthresh.brier.aa);
mean(mapplus.brier.aa);
mean(map.brier.aa);
mean(raxml.brier.aa);
mean(post.brier.aa);
mean(rand.brier.aa);







# Varying tuning parameters

source("/home/guindon/latex/student/adrien/revision/res/ancestral_analyse.R");
filename="/home/guindon/latex/student/adrien/revision/res/merged_2018_12_26/Codon-900-withProba/all";
d=read.table(filename,header=T,comment="#");
p = as.matrix(cbind(d$pA,d$pC,d$pG,d$pT));
s = pbapply(p,1,map.plus,0.05);
d$mapplus=s;
s = pbapply(p,1,cum.prob,0.5);
d$cumprob=s;
s = pbapply(p,1,prob.thresh,0.45);
d$probthresh=s;
d = d[apply(cbind(d$pA,d$pC,d$pG,d$pT),1,max) < 0.95,];
mapplus.accu = accuracy(as.vector(d$True),as.vector(d$mapplus));
mapplus.prec = precision(as.vector(d$mapplus));

cumprob.accu = accuracy(as.vector(d$True),as.vector(d$cumprob));
cumprob.prec = precision(as.vector(d$cumprob));

probthresh.accu = accuracy(as.vector(d$True),as.vector(d$probthresh));
probthresh.prec = precision(as.vector(d$probthresh));

mapplus.accu = accuracy(as.vector(d$True),as.vector(d$mapplus));
mapplus.prec = precision(as.vector(d$mapplus));
 
cumprob.accu = accuracy(as.vector(d$True),as.vector(d$cumprob));
cumprob.prec = precision(as.vector(d$cumprob));
 
probthresh.accu = accuracy(as.vector(d$True),as.vector(d$probthresh));
probthresh.prec = precision(as.vector(d$probthresh));
c11 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d11 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e11 = sum(probthresh.accu[probthresh.prec == amb] == 1);

amb=1;
c11 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d11 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e11 = sum(probthresh.accu[probthresh.prec == amb] == 1);

c01 = sum(mapplus.accu[mapplus.prec == amb] == 0);
d01 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e01 = sum(probthresh.accu[probthresh.prec == amb] == 0);

amb=2;
c12 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d12 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e12 = sum(probthresh.accu[probthresh.prec == amb] == 1);

c02 = sum(mapplus.accu[mapplus.prec == amb] == 0);
d02 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e02 = sum(probthresh.accu[probthresh.prec == amb] == 0);

amb=3;
c13 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d13 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e13 = sum(probthresh.accu[probthresh.prec == amb] == 1);

c03 = sum(mapplus.accu[mapplus.prec == amb] == 0);
d03 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e03 = sum(probthresh.accu[probthresh.prec == amb] == 0);

amb = 4;
c14 = sum(mapplus.accu[mapplus.prec >= amb] == 1);
d14 = sum(cumprob.accu[cumprob.prec >= amb] == 1);
e14 = sum(probthresh.accu[probthresh.prec >= amb] == 1);

c04 = sum(mapplus.accu[mapplus.prec >= amb] == 0);
d04 = sum(cumprob.accu[cumprob.prec >= amb] == 0);
e04 = sum(probthresh.accu[probthresh.prec >= amb] == 0);

cat(c01," ",c11," ",c01/(c11+c01)," ",(c11+c01)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d01," ",d11," ",d01/(d11+d01)," ",(d11+d01)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e01," ",e11," ",e01/(e11+e01)," ",(e11+e01)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c01," ",c11," ",c01/(c11+c01)," ",(c11+c01)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d01," ",d11," ",d01/(d11+d01)," ",(d11+d01)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e01," ",e11," ",e01/(e11+e01)," ",(e11+e01)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c02," ",c12," ",c02/(c12+c02)," ",(c12+c02)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d02," ",d12," ",d02/(d12+d02)," ",(d12+d02)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e02," ",e12," ",e02/(e12+e02)," ",(e12+e02)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c03," ",c13," ",c03/(c13+c03)," ",(c13+c03)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d03," ",d13," ",d03/(d13+d03)," ",(d13+d03)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e03," ",e13," ",e03/(e13+e03)," ",(e13+e03)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c04," ",c14," ",c04/(c14+c04)," ",(c14+c04)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d04," ",d14," ",d04/(d14+d04)," ",(d14+d04)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e04," ",e14," ",e04/(e14+e04)," ",(e14+e04)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");






source("/home/guindon/latex/student/adrien/revision/res/ancestral_analyse.R");
filename="/home/guindon/latex/student/adrien/revision/res/merged_2018_12_26/Codon-900-withProba/all";
d=read.table(filename,header=T,comment="#");
p = as.matrix(cbind(d$pA,d$pC,d$pG,d$pT));
s = pbapply(p,1,map.plus,0.30);
d$mapplus=s;
s = pbapply(p,1,cum.prob,0.80);
d$cumprob=s;
s = pbapply(p,1,prob.thresh,0.35);
d$probthresh=s;
d = d[apply(cbind(d$pA,d$pC,d$pG,d$pT),1,max) < 0.95,];
mapplus.accu = accuracy(as.vector(d$True),as.vector(d$mapplus));
mapplus.prec = precision(as.vector(d$mapplus));

cumprob.accu = accuracy(as.vector(d$True),as.vector(d$cumprob));
cumprob.prec = precision(as.vector(d$cumprob));

probthresh.accu = accuracy(as.vector(d$True),as.vector(d$probthresh));
probthresh.prec = precision(as.vector(d$probthresh));

mapplus.accu = accuracy(as.vector(d$True),as.vector(d$mapplus));
mapplus.prec = precision(as.vector(d$mapplus));
 
cumprob.accu = accuracy(as.vector(d$True),as.vector(d$cumprob));
cumprob.prec = precision(as.vector(d$cumprob));
 
probthresh.accu = accuracy(as.vector(d$True),as.vector(d$probthresh));
probthresh.prec = precision(as.vector(d$probthresh));
c11 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d11 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e11 = sum(probthresh.accu[probthresh.prec == amb] == 1);

amb=1;
c11 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d11 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e11 = sum(probthresh.accu[probthresh.prec == amb] == 1);

c01 = sum(mapplus.accu[mapplus.prec == amb] == 0);
d01 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e01 = sum(probthresh.accu[probthresh.prec == amb] == 0);

amb=2;
c12 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d12 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e12 = sum(probthresh.accu[probthresh.prec == amb] == 1);

c02 = sum(mapplus.accu[mapplus.prec == amb] == 0);
d02 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e02 = sum(probthresh.accu[probthresh.prec == amb] == 0);

amb=3;
c13 = sum(mapplus.accu[mapplus.prec == amb] == 1);
d13 = sum(cumprob.accu[cumprob.prec == amb] == 1);
e13 = sum(probthresh.accu[probthresh.prec == amb] == 1);

c03 = sum(mapplus.accu[mapplus.prec == amb] == 0);
d03 = sum(cumprob.accu[cumprob.prec == amb] == 0);
e03 = sum(probthresh.accu[probthresh.prec == amb] == 0);

amb = 4;
c14 = sum(mapplus.accu[mapplus.prec >= amb] == 1);
d14 = sum(cumprob.accu[cumprob.prec >= amb] == 1);
e14 = sum(probthresh.accu[probthresh.prec >= amb] == 1);

c04 = sum(mapplus.accu[mapplus.prec >= amb] == 0);
d04 = sum(cumprob.accu[cumprob.prec >= amb] == 0);
e04 = sum(probthresh.accu[probthresh.prec >= amb] == 0);

cat(c01," ",c11," ",c01/(c11+c01)," ",(c11+c01)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d01," ",d11," ",d01/(d11+d01)," ",(d11+d01)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e01," ",e11," ",e01/(e11+e01)," ",(e11+e01)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c01," ",c11," ",c01/(c11+c01)," ",(c11+c01)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d01," ",d11," ",d01/(d11+d01)," ",(d11+d01)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e01," ",e11," ",e01/(e11+e01)," ",(e11+e01)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c02," ",c12," ",c02/(c12+c02)," ",(c12+c02)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d02," ",d12," ",d02/(d12+d02)," ",(d12+d02)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e02," ",e12," ",e02/(e12+e02)," ",(e12+e02)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c03," ",c13," ",c03/(c13+c03)," ",(c13+c03)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d03," ",d13," ",d03/(d13+d03)," ",(d13+d03)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e03," ",e13," ",e03/(e13+e03)," ",(e13+e03)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");

cat(c04," ",c14," ",c04/(c14+c04)," ",(c14+c04)/(c11+c01+c12+c02+c13+c03+c14+c04),"\n");
cat(d04," ",d14," ",d04/(d14+d04)," ",(d14+d04)/(d11+d01+d12+d02+d13+d03+d14+d04),"\n");
cat(e04," ",e14," ",e04/(e14+e04)," ",(e14+e04)/(e11+e01+e12+e02+e13+e03+e14+e04),"\n");




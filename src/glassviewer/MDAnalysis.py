import glassviewer as pc
import glassviewer.traj_process as ptp
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import time
import pickle
#from multiprocessing import Pool
from multiprocessing import get_context
printlog=1
def maxcluster(sys1,clusterneimethod,clustercutoff,clusterq6threshold,onlyreturnsolidatomsnum,qaverage):
    sys1.find_neighbors(method=clusterneimethod,cutoff='sann')#voronoi 和 cutoff/scann 之间可选
    sys1.calculate_q([6],averaged=qaverage)
    BOO=sys1.get_qvals([6],averaged=qaverage)
    #print('the BOO number (q2, q4, q5, q6) for Atom 0 is '+str(BOO[0][0])+' '+str(BOO[0][1])+' '+str(BOO[0][2])+' '+str(BOO[0][3]))
    def aq6condition(atom):
        if atom.allaq[6-2]>clusterq6threshold:
            return True
        else:
            return False
    #log(str(onlyreturnsolidatomsnum),printlog)
    if onlyreturnsolidatomsnum:
        count=0
        for a in sys1.atoms:
            if aq6condition(a):
                count+=1
        return count
    else:
        return sys1.cluster_atoms(condition=aq6condition,cutoff=clustercutoff)
def log(a,flag):
    if flag:
        with open('Anarunlog.log','a') as flog:
            flog.write('%-40s'%a+'\t'+time.asctime( time.localtime(time.time()))+'\n')
def mkdir(strs):
    newdir='./'+strs+'_dir'
    if os.name=='nt':
        newdir=newdir.replace('/','\\')
        os.system('rmdir /Q /S '+newdir)
    elif os.name=='posix':
        os.system('rm -r '+newdir)
    if not os.path.exists(newdir):
        os.makedirs(newdir,exist_ok=True)
    elif not os.path.exists(newdir+'/'):
        raise ValueError('file with name '+strs+'_dir exists, folder with the same name cannot be created')
def savedata(strs,i,data):
    f=open('./'+strs+'_dir/'+strs+str(i),'wb')
    pickle.dump(data,f)
    f.close()
def readdata(strs,i):
    f=open('./'+strs+'_dir/'+strs+str(i),'rb')
    data=pickle.load(f)
    f.close()
    return data
class MDAnalysis:

    def __init__(self):
        #self.filesname=[]
        self.filelist=['WT_BLPSXDATCAR','WT_NLPSXDATCAR','XWM_BLPSXDATCAR','XWM_NLPSXDATCAR','VASPXDATCAR']
        self.partial=[[1,1],[1,2],[2,1],[2,2]]
        self.partialname=['LiLi','LiMg','MgLi','MgMg']
        self.structurename= ['others', 'fcc', 'hcp', 'bcc', 'ico']
        self.BOOsname= range(2,13)
        self.VoronoiIndexname= range(3,7)
        self.format='poscar'
        self.taskcount=0
        self.processnum=0
        
        self.stepnumber=1
        self.jumpednumber=0
        self.scale=1
        
        self.neighboron=True
        self.pdfon=True
        self.sfon=True
        self.badon=True
        self.SROon=True
        self.CNAon=True
        self.BOOon=True
        self.VIon=True
        self.savefig=False
        self.clusteron=True
        self.BOOtimeson=True
        self.firstpasson=True
        self.Analyze2don=True

        self.plotextraname=''

        self.neighbormethod='voronoi'
        
        self.pdfBins=100
        self.pdfcut=10
        self.pdfthreadnum=1
        self.pdfr=[]
        self.partialpdfs=[]
        self.totalpdfs=[]

        self.sfpdfBins=100
        self.sfpdfcut=10
        self.sfpdfthreadnum=1
        self.sfr=[]
        self.sftotals=[]
        self.sfpartials=[]
        self.sfview_pointcut=0

        self.badBins=100 
        self.bads=[]
        self.badr=[]
        
        self.SRO_Cowleys=[]
        self.SRO_CS_unnorms=[]
        self.SRO_CS_norms=[]
        self.CNAs=[]
        
        self.BOOtimes_Won=False
        self.BOOBins=100
        self.BOOs=[]
        self.BOOr=[]
        self.BOOMax=[]
        self.BOOMin=[]
        self.BOOaverage=True
        
        self.BOOtimesBins=1000
        self.BOOstimes=[]
        #self.BOOrtimes=[]
        self.BOOtimesMax=[]
        self.BOOtimesMin=[]
        self.BOOtimesaverage=True

        #self.VoronoiIndexBins=100
        self.VoronoiIndex=[]
        self.VIr=[]
        self.VoronoiIndexMax=[]
        self.VoronoiIndexMin=[]
            
        self.clusterneimethod='cutoff'
        self.clustercutoff=0
        self.clusterq6threshold=3.5
        self.cluster=[]
        self.clusterqaverage=True
        self.onlyreturnsolidatomsnum=False
        self.firstpasstime=[]
        self.atomsnum=0
        
        self.bins2d=100
        self.cut2d=[3,4]
        self.histq4w4range=[[0,0.25],[-0.2,0.2]]
        self.histq4q8range=[[0,0.25],[0.2,0.5]]
        self.histq4aq6range=[[0,0.25],[0,0.55]]
        self.histaq4aq6range=[[0,0.2],[0,0.6]]
        self.neighbordistrange=(10,18)
        self.histq4w4=[]
        self.histq4q8=[]
        self.histq4aq6=[]
        self.histaq4aq6=[]
        self.neighbordist=[]
        self.solidnum=0
        self.precursornum=0
        self.neighborthreshold2d=2
        self.precursoraq6cut=0.2
        self.showmode='all' #all precursor solid
    def turnoffall(self):
        self.neighboron=False
        self.pdfon=False
        self.sfon=False
        self.badon=False
        self.SROon=False
        self.CNAon=False
        self.BOOon=False
        self.VIon=False
        self.clusteron=False
        self.BOOtimeson=False
        self.firstpasson=False
        self.Analyze2don=False
    def MDplot(self):
        #sns.set()
        plt.style.use(['science','no-latex'])
        
        if self.Analyze2don:

            #basecolorbar=np.array([0.02,0.05,0.1,0.2,0.3,0.5,0.7,0.9,1])
            #basecolorbar=np.array(np.linspace(0,1,11))
            def colorbarfunc(array):
                array=np.array(array)
                return np.linspace(array[~np.isneginf(array)].min(),array[~np.isneginf(array)].max(),11)
            
            with np.errstate(divide='ignore'):
                self.histq4q8log=np.log(self.histq4q8)
                self.histq4w4log=np.log(self.histq4w4)
                self.histq4aq6log=np.log(self.histq4aq6)
                self.histaq4aq6log=np.log(self.histaq4aq6)

                
            heightbarq4q8=colorbarfunc(self.histq4q8log)
            heightbarq4w4=colorbarfunc(self.histq4w4log)
            heightbarq4aq6=colorbarfunc(self.histq4aq6log)
            heightbaraq4aq6=colorbarfunc(self.histaq4aq6log)
            
            self.histq4q8log[np.isneginf(self.histq4q8log)]=1e100
            self.histq4w4log[np.isneginf(self.histq4w4log)]=1e100
            self.histq4aq6log[np.isneginf(self.histq4aq6log)]=1e100
            self.histaq4aq6log[np.isneginf(self.histaq4aq6log)]=1e100
            
            #from IPython import embed
            #embed()
            hist,X,Y=np.histogram2d(x=[self.histq4q8range[0][0]],y=[self.histq4q8range[1][0]],bins=self.bins2d,density=True,range=self.histq4q8range)
            x,y=np.meshgrid(X[0:len(X)-1],Y[0:len(Y)-1],indexing='ij')
            #heightbar=colorbarfunc(self.histq4q8)
            plt.figure(dpi=150)
            sns.set()
            plt.contourf(x,y,self.histq4q8log,heightbarq4q8,cmap='tab20b')
            plt.colorbar()
            plt.title('q4-q8 distribution')
            plt.xlabel('q4')
            plt.ylabel('q8')
            if(self.savefig==True):
                plt.savefig('q4-q8'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
            
            hist,X,Y=np.histogram2d(x=[self.histq4w4range[0][0]],y=[self.histq4w4range[1][0]],bins=self.bins2d,density=True,range=self.histq4w4range)
            x,y=np.meshgrid(X[0:len(X)-1],Y[0:len(Y)-1],indexing='ij')
            #heightbar=colorbarfunc(self.histq4w4)
            plt.figure(dpi=150)
            sns.set()
            #from IPython import embed
            #embed()
            plt.contourf(x,y,self.histq4w4log,heightbarq4w4,cmap='tab20b')
            plt.colorbar()
            plt.title('q4-w4 distribution')
            plt.xlabel('q4')
            plt.ylabel('w4')
            if(self.savefig==True):
                plt.savefig('q4-w4'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
            
            hist,X,Y=np.histogram2d(x=[self.histq4aq6range[0][0]],y=[self.histq4aq6range[1][0]],bins=self.bins2d,density=True,range=self.histq4aq6range)
            x,y=np.meshgrid(X[0:len(X)-1],Y[0:len(Y)-1],indexing='ij')
            #heightbar=colorbarfunc(self.histq4aq6)
            plt.figure(dpi=150)
            sns.set()
            plt.contourf(x,y,self.histq4aq6log,heightbarq4aq6,cmap='tab20b')
            plt.colorbar()
            plt.title('q4-aq6 distribution')
            plt.xlabel('q4')
            plt.ylabel('aq6')
            if(self.savefig==True):
                plt.savefig('q4-aq6'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
            
            hist,X,Y=np.histogram2d(x=[self.histaq4aq6range[0][0]],y=[self.histaq4aq6range[1][0]],bins=self.bins2d,density=True,range=self.histaq4aq6range)
            x,y=np.meshgrid(X[0:len(X)-1],Y[0:len(Y)-1],indexing='ij')
            #heightbar=colorbarfunc(self.histaq4aq6)
            plt.figure(dpi=150)
            sns.set()
            plt.contourf(x,y,self.histaq4aq6log,heightbaraq4aq6,cmap='tab20b')
            plt.colorbar()
            plt.title('aq4-aq6 distribution')
            plt.xlabel('aq4')
            plt.ylabel('aq6')
            if(self.savefig==True):
                plt.savefig('aq4-aq6'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
    
            plt.figure(dpi=150)
            plt.bar(x=range(self.neighbordistrange[0],self.neighbordistrange[1]),height=self.neighbordist)
            plt.title('neighbor number distribution')
            plt.xlabel('neighbor number')
            plt.ylabel('count')
            if(self.savefig==True):
                plt.savefig('neighbornumberdistribution'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
            
            
        if self.clusteron:
            for i,x in enumerate(self.filelist):
                plt.figure(dpi=150)
                sns.set()
                plt.title(self.filelist[i] +' ClusterSize')
                plt.xlabel('x')
                plt.ylabel('y')
                plt.plot(range(1,self.stepnumber+1),self.cluster[i][::-1])#由于倒序读的，所以调换顺序
                if(self.savefig==True):
                    plt.savefig(self.filelist[i] + 'ClusterSize'+self.plotextraname+'.png')
                    plt.close()
                else:
                    plt.show()
        if self.firstpasson:
            for i,x in enumerate(self.filelist):
                plt.figure(dpi=150)
                sns.set()
                plt.title(self.filelist[i] +' FirstPassTime')
                plt.xlabel('x')
                plt.ylabel('y')
                plt.plot(range(1,self.atomsnum+1),self.firstpasstime[i])
                if(self.savefig==True):
                    plt.savefig(self.filelist[i] + 'FirstPassTime'+self.plotextraname+'.png')
                    plt.close()
                else:
                    plt.show()
        if self.pdfon:
            for i,x in enumerate(self.partial):
                plt.figure(dpi=150)
                sns.set()
                plt.title(self.partialname[i] +' partialPDF')
                plt.xlabel('x')
                plt.ylabel('y')
                for j in range(len(self.filelist)):
                    plt.plot(self.pdfr,self.partialpdfs[i][j])
                plt.legend(self.filelist)
                if(self.savefig==True):
                    plt.savefig(self.partialname[i] + 'partialPDF'+self.plotextraname+'.png')
                    plt.close()
                else:
                    plt.show()

            plt.figure(dpi=150)
            sns.set()
            plt.title('TotalPDF')
            plt.xlabel('x')
            plt.ylabel('y')
            for i in range(len(self.filelist)):
                plt.plot(self.pdfr,self.totalpdfs[i])
            plt.legend(self.filelist)
            if(self.savefig==True):
                plt.savefig('TotalPDF'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
        if self.sfon:
            for i,x in enumerate(self.partial):
                plt.figure(dpi=150)
                sns.set()
                plt.title(self.partialname[i] +' partialSF')
                plt.xlabel('x')
                plt.ylabel('y')
                for j in range(len(self.filelist)):
                    if self.sfview_pointcut==0:
                        plt.plot(self.sfr,self.sfpartials[i][j])
                    else:
                        plt.plot(self.sfr[0:self.sfview_pointcut],self.sfpartials[i][j][0:self.sfview_pointcut])
                plt.legend(self.filelist)
                if(self.savefig==True):
                    plt.savefig(self.partialname[i] +' partialSF'+self.plotextraname+'.png')
                    plt.close()
                else:
                    plt.show()

            plt.figure(dpi=150)
            sns.set()
            plt.title('TotalSF')
            plt.xlabel('x')
            plt.ylabel('y')
            for i in range(len(self.filelist)):
                if self.sfview_pointcut==0:
                    plt.plot(self.sfr,self.sftotals[i])
                else:
                    plt.plot(self.sfr[0:self.sfview_pointcut],self.sftotals[i][0:self.sfview_pointcut])
            plt.legend(self.filelist)
            if(self.savefig==True):
                plt.savefig('TotalSF'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
        if self.badon:
            plt.figure(dpi=150)
            sns.set()
            plt.title('BAD')
            plt.xlabel('x')
            plt.ylabel('y')
            for i in range(len(self.filelist)):
                plt.plot(self.badr,self.bads[i])
            plt.legend(self.filelist)
            if(self.savefig==True):
                plt.savefig('BAD'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
        if self.SROon:
            plt.figure(dpi=150)
            sns.set()
            for i,x in enumerate(self.partial):
                plt.scatter(np.arange(len(self.SRO_Cowleys[i])),self.SRO_Cowleys[i])
            plt.xticks(range(len(self.filelist)),self.filelist, rotation=30, fontsize=10)
            plt.title('SRO_Cowley')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.legend(self.partialname)
            if(self.savefig==True):
                plt.savefig('SRO_Cowley'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()

            plt.figure(dpi=150)
            sns.set()
            for i,x in enumerate(self.partial):
                plt.scatter(np.arange(len(self.SRO_CS_unnorms[i])),self.SRO_CS_unnorms[i])
            plt.xticks(range(len(self.filelist)),self.filelist, rotation=30, fontsize=10)
            plt.title('SRO_CS_unnorms')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.legend(self.partialname)
            if(self.savefig==True):
                plt.savefig('SRO_CS_unnorms'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()

            plt.figure(dpi=150)
            sns.set()
            for i,x in enumerate(self.partial):
                plt.scatter(np.arange(len(self.SRO_CS_norms[i])),self.SRO_CS_norms[i])
            plt.xticks(range(len(self.filelist)),self.filelist, rotation=30, fontsize=10)
            plt.title('SRO_CS_norms')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.legend(self.partialname)
            if(self.savefig==True):
                plt.savefig('SRO_CS_norms'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
        if self.CNAon:
            plt.figure(dpi=150)
            sns.set()
            for i,x in enumerate(self.filelist):
                plt.plot(self.CNAs[i])
            plt.xticks(range(len(self.structurename)),self.structurename, rotation=30, fontsize=10)
            plt.title('CNA')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.legend(self.filelist)
            if(self.savefig==True):
                plt.savefig('CNA'+self.plotextraname+'.png')
                plt.close()
            else:
                plt.show()
        if self.BOOon:
            for p,x in enumerate(self.BOOsname):
                plt.figure(dpi=150)
                sns.set()
                plt.title('BOO_Q'+str(x))
                plt.xlabel('x')
                plt.ylabel('y')
                for i in range(len(self.filelist)):
                    plt.plot(self.BOOr[p],self.BOOs[p][i])
                plt.legend(self.filelist)
                if(self.savefig==True):
                    plt.savefig('BOO_Q'+str(x)+''+self.plotextraname+'.png')
                    plt.close()
                else:
                    plt.show()
        if self.BOOtimeson:
            
            def colorbarfunc(array):
                array=np.array(array)
                return np.linspace(array[~np.isneginf(array)].min(),array[~np.isneginf(array)].max(),11)

            BOOstimestemp2=np.array([[[self.BOOstimes[-1-i][Pid][Xid] for i in range(self.stepnumber)] for Xid in range(len(self.filelist))] for Pid in range(len(self.BOOsname))])
            for ifile,xfile in enumerate(self.filelist):
                for iboo,xboo in enumerate(self.BOOsname):
                    f,ax=plt.subplots(figsize=(5,5/180*self.stepnumber),dpi=200)
                    
                    x=np.linspace(self.BOOtimesMin[iboo],self.BOOtimesMax[iboo],self.BOOtimesBins,endpoint=False)
                    y=np.arange(0,self.stepnumber)
                    X,Y=np.meshgrid(x,y)
                    
                    out=BOOstimestemp2[iboo][ifile]
                    with np.errstate(divide='ignore'):
                        out=np.log(out)       
                    heightbar=colorbarfunc(out)
                    out[np.isneginf(out)]=1e100
                    
                    #heightbar=np.array([0.02,0.05,0.1,0.2,0.3,0.5,0.7,0.9,1])*out.max()
                    m=ax.contourf(X,Y,out,heightbar,cmap='tab20b')#
                    ax.contour(X,Y,out,heightbar,colors='black',linewidths=0.5)
                    
                    plt.colorbar(m,ax=ax)
                    
                    #plt.minorticks_on()
                    if not self.BOOtimes_Won:
                        if self.BOOtimesaverage:
                            plt.title('Boo Average Q'+str(xboo)+' Value Distribution')
                            plt.xlabel('average q'+str(xboo))
                        else:
                            plt.title('Boo Q'+str(xboo)+' Value Distribution')
                            plt.xlabel('q'+str(xboo))
                        plt.ylabel('Frame')
                        if(self.savefig==True):
                            if self.BOOtimesaverage:
                                plt.savefig('boo_vs_time_aq'+str(xboo)+'_'+xfile+self.plotextraname+'.png')
                            else:
                                plt.savefig('boo_vs_time_q'+str(xboo)+'_'+xfile+self.plotextraname+'.png')
                        else:
                            plt.show()
                    else:
                        if self.BOOtimesaverage:
                            plt.title('Boo Average W_norm'+str(xboo)+' Value Distribution')
                            plt.xlabel('average W_norm'+str(xboo))
                        else:
                            plt.title('Boo W_norm'+str(xboo)+' Value Distribution')
                            plt.xlabel('W_norm'+str(xboo))
                        plt.ylabel('Frame')
                        if(self.savefig==True):
                            if self.BOOtimesaverage:
                                plt.savefig('boo_vs_time_awnorm'+str(xboo)+'_'+xfile+self.plotextraname+'.png')
                            else:
                                plt.savefig('boo_vs_time_worm'+str(xboo)+'_'+xfile+self.plotextraname+'.png')
                        else:
                            plt.show()
        if self.VIon:
            for p,x in enumerate(self.VoronoiIndexname):
                plt.figure(dpi=150)
                sns.set()
                plt.title('Voronoi Index'+str(x))
                plt.xlabel('Voronoi index'+str(x)+' Value')
                plt.ylabel('Averaged count')
                width=1.0/(len(self.filelist)+1)
                for i in range(len(self.filelist)):
                    plt.bar(np.array(self.VIr[p])+width*(i-(len(self.filelist)-1)/2.0),self.VoronoiIndex[p][i],width=width)
                plt.xticks(self.VIr[p])
                plt.legend(self.filelist,bbox_to_anchor=(1,1))
                if(self.savefig==True):
                    plt.savefig('Voronoi Index'+str(x)+''+self.plotextraname+'.png')
                    plt.close()
                else:
                    plt.show()
    def getnecepara(self,filesname):
        
        for XDATCARNo in range(len(self.filelist)):
            if XDATCARNo==0:
                sys=pc.System()
                sys.read_inputfile(filesname[XDATCARNo][0-self.jumpednumber-1],self.format)
                if self.neighboron:
                    if self.neighbormethod=='voronoi':
                        sys.find_neighbors(method=self.neighbormethod)
                    elif self.neighbormethod=='cutoff-sann':
                        sys.find_neighbors(method='cutoff',cutoff='sann')
                    elif self.neighbormethod=='cutoff-fit':
                        import warnings
                        warnings.simplefilter('ignore', np.RankWarning)
                        pdfbintemp=int(200*(self.cut2d[1]-self.cut2d[0]))
                        
                        pdf, r=sys.calculate_pdf(histobins=pdfbintemp,histomin=self.cut2d[0],cut=self.cut2d[1],threadnum=1)
                        c1=np.polyfit(r,pdf,20)
                        yfit1=np.polyval(c1,r)
                        cut1=(np.argmin(yfit1))
                        
                        min2=np.max([0,cut1-100])
                        max2=np.min([pdfbintemp,cut1+100])
                        
                        c2=np.polyfit(r[min2:max2],pdf[min2:max2],20)
                        yfit2=np.polyval(c2,r)
                        cut2=(np.argmin(yfit2[min2:max2])+min2)
                        cutoffv=self.cut2d[0]+cut2*(self.cut2d[1]-self.cut2d[0])/pdfbintemp
                        
                        sys.find_neighbors(method='cutoff',cutoff=cutoffv)
                for j,x in enumerate(self.partial):
                    if self.pdfon:
                        partialpdf, self.pdfr=sys.calculate_pdf(histobins=self.pdfBins,histomin=0,cut=self.pdfcut,threadnum=self.pdfthreadnum,partial=True,centertype=x[0],secondtype=x[1])
                        #self.partialpdfs[j][XDATCARNo]=partialpdf
                    if self.sfon:
                        sfpartialpdf, sfpdfr=sys.calculate_pdf(histobins=self.sfpdfBins,histomin=0,cut=self.sfpdfcut,threadnum=self.sfpdfthreadnum,partial=True,centertype=x[0],secondtype=x[1])
                        sfpartial,self.sfr=sys.calculate_sf(sfpartialpdf, sfpdfr,0);
                        #self.sfpartials[j][XDATCARNo]=np.array(sfpartial)
                   # if self.SROon:
                        #self.SRO_Cowleys[j][XDATCARNo]=sys.calculate_pmsro(reference_type=x[0],compare_type=x[1])[0]
                        #self.SRO_CS_unnorms[j][XDATCARNo]=sys.calculate_pmsro_CS(reference_type=x[0],compare_type=x[1],normalization=False)
                        #self.SRO_CS_norms[j][XDATCARNo]=sys.calculate_pmsro_CS(reference_type=x[0],compare_type=x[1],normalization=True)
                if self.pdfon:
                    totalpdf, self.pdfr=sys.calculate_pdf(histobins=self.pdfBins,histomin=0,cut=self.pdfcut,threadnum=self.pdfthreadnum)
                    #self.totalpdfs[XDATCARNo]=totalpdf
                if self.sfon:
                    sftotalpdf, sfpdfr=sys.calculate_pdf(histobins=self.sfpdfBins,histomin=0,cut=self.sfpdfcut,threadnum=self.sfpdfthreadnum)
                    sftotal,self.sfr=sys.calculate_sf(sftotalpdf, sfpdfr,0);
                    #self.sftotals[XDATCARNo]=np.array(sftotal)
                if self.badon:
                    bad, self.badr=sys.calculate_bad(histobins=self.badBins,histomin=0,histomax=np.pi);
                    #self.bads[XDATCARNo]=bad
                #if self.CNAon:
                    #self.CNAs[XDATCARNo]=list(sys.calculate_cna().values())
                if self.BOOon:
                    sys.calculate_q(self.BOOsname,averaged=self.BOOaverage)
                    BOO=sys.get_qvals(self.BOOsname,averaged=self.BOOaverage)#2-13
                    for p,x in enumerate(self.BOOsname):
                        self.BOOMax[p]=max(BOO[p])
                        self.BOOMin[p]=min(BOO[p])
                        BOOhist,BOOrtemp=np.histogram(BOO[p],range=(self.BOOMin[p],self.BOOMax[p]),bins=self.BOOBins)
                        #self.BOOs[p][XDATCARNo]=BOOhist
                        self.BOOr[p]=BOOrtemp[0:self.BOOBins]
                if self.BOOtimeson:
                    pass
                if self.VIon:
                    sys.calculate_vorovector()#3 4 5 6
                    for p,x in enumerate(self.VoronoiIndexname):
                        VI=[a.vorovector[p] for a in sys.atoms]
                        self.VoronoiIndexMax[p]=int(max(VI))
                        self.VoronoiIndexMin[p]=int(min(VI))
                        for tp in range(len(self.filelist)):
                            self.VoronoiIndex[p][tp]=np.zeros(int(self.VoronoiIndexMax[p]+1-self.VoronoiIndexMin[p]))
                        VIhist,VIrtemp=np.histogram(VI,range=(self.VoronoiIndexMin[p]-0.5,self.VoronoiIndexMax[p]+0.5),bins=self.VoronoiIndexMax[p]+1-self.VoronoiIndexMin[p])
                        #self.VoronoiIndex[p][XDATCARNo]=VIhist
                        self.VIr[p]=range(self.VoronoiIndexMin[p],self.VoronoiIndexMax[p]+1)
                if self.firstpasson:
                    self.atomsnum=sys.natoms
    def callback(self,result):
        self.taskcount+=1
        log('%5d/%5d\tPer%10.2f'%(self.taskcount,self.stepnumber,self.taskcount/self.stepnumber*100)+'%',printlog)

def calculator(MD,filesname):
        #初始化 
    pool=get_context('spawn').Pool(processes=MD.processnum)
    #pool=Pool(processes=MD.processnum)
    if MD.Analyze2don:
        MD.histq4w4=np.zeros((MD.bins2d,MD.bins2d))
        MD.histq4q8=np.zeros((MD.bins2d,MD.bins2d))
        MD.histq4aq6=np.zeros((MD.bins2d,MD.bins2d))
        MD.histaq4aq6=np.zeros((MD.bins2d,MD.bins2d))
        MD.neighbordist=np.zeros(MD.neighbordistrange[1]-MD.neighbordistrange[0])
        mkdir('analyze2d')
    if MD.clusteron:
        MD.cluster=np.array([np.zeros(MD.stepnumber)for i in range(len(MD.filelist))])
        mkdir('cluster')
    if MD.pdfon:
        MD.pdfr=np.zeros(MD.pdfBins) 
        MD.partialpdfs=[[np.zeros(MD.pdfBins) for i in range(len(MD.filelist))]for i in MD.partial]
        MD.totalpdfs=[np.zeros(MD.pdfBins) for i in range(len(MD.filelist))]
        mkdir('pdf')
    if MD.sfon:
        MD.sfr=np.zeros(MD.sfpdfBins) 
        MD.sfpartials=[[np.zeros(MD.sfpdfBins) for i in range(len(MD.filelist))]for i in MD.partial]
        MD.sftotals=[np.zeros(MD.sfpdfBins) for i in range(len(MD.filelist))]
        mkdir('sf')
    if MD.badon:
        MD.bads=[np.zeros(MD.badBins) for i in range(len(MD.filelist))]
        mkdir('bad')
    if MD.SROon:
        MD.SRO_Cowleys=[[0 for i in range(len(MD.filelist))]for i in MD.partial]
        MD.SRO_CS_unnorms=[[0 for i in range(len(MD.filelist))]for i in MD.partial]
        MD.SRO_CS_norms=[[0 for i in range(len(MD.filelist))]for i in MD.partial]
        mkdir('sro')
    if MD.CNAon:
        MD.CNAs=[np.zeros(len(MD.structurename)) for i in range(len(MD.filelist))]
        mkdir('cna')
    if MD.BOOon:
        MD.BOOs=[[np.zeros(MD.BOOBins) for i in range(len(MD.filelist))] for i in MD.BOOsname]
        MD.BOOr=[np.zeros(MD.BOOBins) for i in MD.BOOsname]
        MD.BOOMax=[0  for i in MD.BOOsname]
        MD.BOOMin=[0  for i in MD.BOOsname]
        mkdir('boo')        
    if MD.BOOtimeson:
        MD.BOOstimestemp=[[np.zeros(MD.BOOtimesBins) for j in range(len(MD.filelist))] for i in MD.BOOsname]
        MD.BOOstimes=[[[np.zeros(MD.BOOtimesBins) for i in range(len(MD.filelist))] for i in MD.BOOsname] for a in range(MD.stepnumber)]
        if MD.BOOtimesMax==[]:
            MD.BOOtimesMax=[0  for i in MD.BOOsname]
            for p,x in enumerate(MD.BOOsname):
                MD.BOOtimesMax[p]=0.7

        if MD.BOOtimesMin==[]:
            MD.BOOtimesMin=[0  for i in MD.BOOsname]
            for p,x in enumerate(MD.BOOsname):
                MD.BOOtimesMin[p]=0
            #BOOhist,BOOrtemp=np.histogram(np.zeros[5],range=(MD.BOOtimesMin[p],MD.BOOtimesMax[p]),bins=MD.BOOtimesBins)
            #MD.BOOrtimes[p]=BOOrtemp[0:MD.BOOtimesBins]
        mkdir('bootimes')
    if MD.VIon:
        MD.VoronoiIndex=[[[] for i in range(len(MD.filelist))] for i in MD.VoronoiIndexname]
        MD.VIr=[[]for i in MD.VoronoiIndexname]
        MD.VoronoiIndexMax=[0  for i in MD.VoronoiIndexname]
        MD.VoronoiIndexMin=[0  for i in MD.VoronoiIndexname]
        mkdir('vi')
    #获得必要的参数
    log('Aquiring necessary parameter for this task',printlog)
    MD.getnecepara(filesname)

    log('MainMPI:begin to distribute MPI tasks',printlog)
    #初始化计数变量
    MD.taskcount=0
    #运行MPI
    #ftemp=open('MDtemp','wb')
    #pickle.dump(MD,ftemp)
    #ftemp.close()
    paralist=[[MD,i,[filesname[XDATCARNo][0-MD.jumpednumber-1-i*MD.scale] for XDATCARNo,x in enumerate(MD.filelist)] ] for i in range(MD.stepnumber)]
    for x in paralist:
        #x=1
        pool.apply_async(calculate_thread,(x,),callback=MD.callback)
    pool.close()
    pool.join()
    
    #对结果求平均
    for i in range(MD.stepnumber):
        if MD.Analyze2don:
            temp=readdata('analyze2d',i)
            MD.solidnum+=temp[0]/MD.stepnumber
            MD.precursornum+=temp[1]/MD.stepnumber
            MD.histq4w4+=temp[2]/MD.stepnumber
            MD.histq4q8+=temp[3]/MD.stepnumber
            MD.histq4aq6+=temp[4]/MD.stepnumber
            MD.histaq4aq6+=temp[5]/MD.stepnumber
            MD.neighbordist+=temp[6]/MD.stepnumber
        if MD.clusteron:
            temp=readdata('cluster',i)
            for a,id in enumerate(temp):
                MD.cluster[a][i]=temp[a]
        if MD.pdfon:
            temp=readdata('pdf',i)
            MD.totalpdfs+=temp[0]/MD.stepnumber
            MD.partialpdfs+=temp[1]/MD.stepnumber
        if MD.sfon:
            temp=readdata('sf',i)
            MD.sftotals+=temp[0]/MD.stepnumber
            MD.sfpartials+=temp[1]/MD.stepnumber
        if MD.badon:
            temp=readdata('bad',i)
            MD.bads+=temp/MD.stepnumber
        if MD.SROon:
            temp=readdata('sro',i)
            MD.SRO_Cowleys+=temp[0]/MD.stepnumber
            MD.SRO_CS_unnorms+=temp[1]/MD.stepnumber
            MD.SRO_CS_norms+=temp[2]/MD.stepnumber
        if MD.CNAon:
            temp=readdata('cna',i)
            MD.CNAs+=temp/MD.stepnumber
        if MD.BOOon:
            temp=readdata('boo',i)
            MD.BOOs+=temp/MD.stepnumber
        if MD.BOOtimeson:
            temp=readdata('bootimes',i)
            MD.BOOstimes[i]=temp
        if MD.VIon:
            temp=readdata('vi',i)
            for a,i in enumerate(temp):
                for b,j in enumerate(i):
                    #print(MD.VoronoiIndex[a][b])
                    #print('#####')
                    #print(j)
                    MD.VoronoiIndex[a][b]+=j/MD.stepnumber
    if MD.firstpasson:
        MD.firstpasstime=[np.zeros(MD.atomsnum) for a in range(len(MD.filelist))]
        for a in range(len(MD.filelist)):
            for b in range(MD.atomsnum):
                for i in range(MD.stepnumber)[::-1]:
                    if(MD.cluster[a][i]>=b+1):
                        MD.firstpasstime[a][b]=MD.stepnumber-i
                        break
    
def calculate_thread(para):
    MD=para[0]
    i=para[1]
    
    filename=para[2]   
        
    

    #for XDATCARNo in range(4):
       # pass
    for XDATCARNo in range(len(MD.filelist)):
        sys=pc.System()
        sys.read_inputfile(filename[XDATCARNo],MD.format)
        if MD.neighboron:
            if MD.neighbormethod=='voronoi':
                sys.find_neighbors(method=MD.neighbormethod)
            elif MD.neighbormethod=='cutoff-sann':
                sys.find_neighbors(method='cutoff',cutoff='sann')
            elif MD.neighbormethod=='cutoff-fit':
                import warnings
                warnings.simplefilter('ignore', np.RankWarning)
                pdfbintemp=int(200*(MD.cut2d[1]-MD.cut2d[0]))
                
                pdf, r=sys.calculate_pdf(histobins=pdfbintemp,histomin=MD.cut2d[0],cut=MD.cut2d[1],threadnum=1)
                c1=np.polyfit(r,pdf,20)
                yfit1=np.polyval(c1,r)
                cut1=(np.argmin(yfit1))
                
                min2=np.max([0,cut1-100])
                max2=np.min([pdfbintemp,cut1+100])
                
                c2=np.polyfit(r[min2:max2],pdf[min2:max2],20)
                yfit2=np.polyval(c2,r)
                cut2=(np.argmin(yfit2[min2:max2])+min2)
                cutoffv=MD.cut2d[0]+cut2*(MD.cut2d[1]-MD.cut2d[0])/pdfbintemp
                
                sys.find_neighbors(method='cutoff',cutoff=cutoffv)
        if MD.Analyze2don:
            
            atom13list=set()
            for atomi,atom in enumerate(sys.atoms):
                if len(atom.neighbors)==13:
                    atom13list.add(atomi)
            
            atom13to12list=set()
            atom13to14list=set()
            if len(atom13list)>0:
                sys.find_neighbors_bynumber_atomlist(threshold=MD.neighborthreshold2d,nmax=14,atomlist=list(atom13list))
            
            atomstemp=sys.atoms
            
            for atomi in atom13list:
                atom=atomstemp[atomi]
                delta= (atom.neighbor_distance[13-1]+atom.neighbor_distance[14-1])/2-cutoffv
                if delta>0:
                    atom13to12list.add(atomi)
                elif delta<=0:
                    atom13to14list.add(atomi)
            if len(atom13to12list)>0:
                sys.find_neighbors_bynumber_atomlist(threshold=MD.neighborthreshold2d,nmax=12,atomlist=list(atom13to12list))
            if len(atom13to14list)>0:
                sys.find_neighbors_bynumber_atomlist(threshold=MD.neighborthreshold2d,nmax=14,atomlist=list(atom13to14list))
            
            sys.calculate_q([4,6,8],averaged=False,clear_condition=True)
            sys.find_solids(bonds=7,threshold=0.7,cluster=False)
            sys.calculate_q([4,6,8],averaged=True,condition='solid',only_averaged=True)
            sys.calculate_w(4,averaged=False)
            
            atoms=set(sys.atoms)
            solidlist=set()
            for atomi,atom in enumerate(atoms):
                if atom.solid:
                    solidlist.add(atom)
            
            liquidlist=set()
            for atomi,atom in enumerate(atoms):
                if not atom.solid:
                    liquidlist.add(atom)
            
            aq6_cutlist=set()
            for atomi,atom in enumerate(atoms):
                if atom.get_q(6,True)>MD.precursoraq6cut:
                    aq6_cutlist.add(atom)
            
            precursorlist=liquidlist.intersection(aq6_cutlist)
            
            if MD.showmode=='precursor':
                showlist=precursorlist
            if MD.showmode=='solid':
                showlist=solidlist
            if MD.showmode=='all':
                showlist=atoms
            
            q4=[]
            q6=[]
            q8=[]
            aq4=[]
            aq6=[]
            w4=[]
            for atom in showlist:
                q4.append(atom.get_q(4,False))
                q6.append(atom.get_q(6,False))
                q8.append(atom.get_q(8,False))
                aq4.append(atom.get_q(4,True))
                aq6.append(atom.get_q(6,True))
                #w4.append(sys.w[0][atom.loc])   
                w4.append(atom.wnorm[4-2])
            MD.solidnum=len(solidlist)
            MD.precursornum=len(precursorlist)
            MD.histq4w4,X,Y=np.histogram2d(x=q4,y=w4,bins=MD.bins2d,density=False,range=MD.histq4w4range)
            MD.histq4q8,X,Y=np.histogram2d(x=q4,y=q8,bins=MD.bins2d,density=False,range=MD.histq4q8range)
            MD.histq4aq6,X,Y=np.histogram2d(x=q4,y=aq6,bins=MD.bins2d,density=False,range=MD.histq4aq6range)
            MD.histaq4aq6,X,Y=np.histogram2d(x=aq4,y=aq6,bins=MD.bins2d,density=False,range=MD.histaq4aq6range)
            a=np.array([len(atom.neighbors) for atom in atoms])
            MD.neighbordist,r=np.histogram(range=MD.neighbordistrange,a=a,bins=(MD.neighbordistrange[1]-MD.neighbordistrange[0]))
            
        if MD.clusteron:
            clustertemp=np.zeros(len(MD.filelist))
            clustertemp[XDATCARNo]=maxcluster(sys,clusterneimethod=MD.clusterneimethod,clustercutoff=MD.clustercutoff,clusterq6threshold=MD.clusterq6threshold,onlyreturnsolidatomsnum=MD.onlyreturnsolidatomsnum,qaverage=MD.clusterqaverage)

        for j,x in enumerate(MD.partial):
            if MD.pdfon:
                partialpdf, pdfr=sys.calculate_pdf(histobins=MD.pdfBins,histomin=0,cut=MD.pdfcut,threadnum=MD.pdfthreadnum,partial=True,centertype=x[0],secondtype=x[1])
                MD.partialpdfs[j][XDATCARNo]=partialpdf

            if MD.sfon:
                sfpartialpdf, sfpdfr=sys.calculate_pdf(histobins=MD.sfpdfBins,histomin=0,cut=MD.sfpdfcut,threadnum=MD.sfpdfthreadnum,partial=True,centertype=x[0],secondtype=x[1])
                sfpartial,sfr=sys.calculate_sf(sfpartialpdf, sfpdfr,0);
                MD.sfpartials[j][XDATCARNo]=np.array(sfpartial)

            if MD.SROon:
                MD.SRO_Cowleys[j][XDATCARNo]=sys.calculate_pmsro(reference_type=x[0],compare_type=x[1])[0]
                MD.SRO_CS_unnorms[j][XDATCARNo]=sys.calculate_pmsro_CS(reference_type=x[0],compare_type=x[1],normalization=False)
                MD.SRO_CS_norms[j][XDATCARNo]=sys.calculate_pmsro_CS(reference_type=x[0],compare_type=x[1],normalization=True)
        if MD.pdfon:
            totalpdf, pdfr=sys.calculate_pdf(histobins=MD.pdfBins,histomin=0,cut=MD.pdfcut,threadnum=MD.pdfthreadnum)
            MD.totalpdfs[XDATCARNo]=totalpdf

        if MD.sfon:
            sftotalpdf, sfpdfr=sys.calculate_pdf(histobins=MD.sfpdfBins,histomin=0,cut=MD.sfpdfcut,threadnum=MD.sfpdfthreadnum)
            sftotal,sfr=sys.calculate_sf(sftotalpdf, sfpdfr,0);
            MD.sftotals[XDATCARNo]=np.array(sftotal)

        if MD.badon:
            bad, badr=sys.calculate_bad(histobins=MD.badBins,histomin=0,histomax=np.pi);
            MD.bads[XDATCARNo]=bad

        if MD.CNAon:
            MD.CNAs[XDATCARNo]=list(sys.calculate_cna().values())
        if MD.BOOon:
            sys.calculate_q(MD.BOOsname,averaged=MD.BOOaverage)
            BOO=sys.get_qvals(MD.BOOsname,averaged=MD.BOOaverage)#2-13
            for p,x in enumerate(MD.BOOsname):
                BOOhist,BOOrtemp=np.histogram(BOO[p],range=(MD.BOOMin[p],MD.BOOMax[p]),bins=MD.BOOBins)
                delta=float(MD.BOOMax[p]-MD.BOOMin[p])/float(MD.BOOBins)
                distri=BOOhist/float(delta*sys.natoms)
                MD.BOOs[p][XDATCARNo]=distri
                #MD.BOOr[p]=BOOrtemp[0:MD.BOOBins]
        if MD.BOOtimeson:
            if not MD.BOOtimes_Won:
                sys.calculate_q(MD.BOOsname,averaged=MD.BOOtimesaverage)
                BOO=sys.get_qvals(MD.BOOsname,averaged=MD.BOOtimesaverage)#2-13
            else:
                sys.calculate_q(MD.BOOsname,averaged=MD.BOOtimesaverage)
                sys.calculate_w(MD.BOOsname,averaged=MD.BOOtimesaverage)
                BOO=[]
                for l in MD.BOOsname:
                    if not MD.BOOtimesaverage:
                        BOO.append([atom.wnorm[l-2] for atom in sys.atoms])
                    else:
                        BOO.append([atom.awnorm[l-2] for atom in sys.atoms])
                        
            for p,x in enumerate(MD.BOOsname):
                BOOhist,BOOrtemp=np.histogram(BOO[p],range=(MD.BOOtimesMin[p],MD.BOOtimesMax[p]),bins=MD.BOOtimesBins)
                delta=float(MD.BOOtimesMax[p]-MD.BOOtimesMin[p])/float(MD.BOOtimesBins)
                distri=BOOhist/float(delta*sys.natoms)
                
                MD.BOOstimestemp[p][XDATCARNo]=distri
                
                #MD.BOOr[p]=BOOrtemp[0:MD.BOOBins]
            
        if MD.VIon:
            sys.calculate_vorovector()#3 4 5 6
            for p,x in enumerate(MD.VoronoiIndexname):
                VI=[a.vorovector[p] for a in sys.atoms]
                VIhist,VIrtemp=np.histogram(VI,range=(MD.VoronoiIndexMin[p]-0.5,MD.VoronoiIndexMax[p]+0.5),bins=MD.VoronoiIndexMax[p]+1-MD.VoronoiIndexMin[p])
                MD.VoronoiIndex[p][XDATCARNo]=VIhist
                #MD.VIr[p]=range(MD.VoronoiIndexMin[p],MD.VoronoiIndexMax[p]+1)
    if MD.Analyze2don:
        savedata('analyze2d',i, [MD.solidnum,MD.precursornum,MD.histq4w4,MD.histq4q8,MD.histq4aq6,MD.histaq4aq6,MD.neighbordist])

    if MD.clusteron:
        savedata('cluster',i, clustertemp)
    if MD.pdfon:
        MD.totalpdfs=np.array(MD.totalpdfs)
        MD.partialpdfs=np.array(MD.partialpdfs)
        savedata('pdf',i, [MD.totalpdfs,MD.partialpdfs])
    if MD.sfon:
        MD.sftotals=np.array(MD.sftotals)
        MD.sfpartials=np.array(MD.sfpartials)
        savedata('sf',i, [MD.sftotals,MD.sfpartials])
    if MD.badon:
        MD.bads=np.array(MD.bads)
        savedata('bad',i, MD.bads)
    if MD.SROon:
        MD.SRO_Cowleys=np.array(MD.SRO_Cowleys)
        MD.SRO_CS_unnorms=np.array(MD.SRO_CS_unnorms)
        MD.SRO_CS_norms=np.array(MD.SRO_CS_norms)
        savedata('sro',i, [MD.SRO_Cowleys,MD.SRO_CS_unnorms,MD.SRO_CS_norms])
    if MD.CNAon:
        MD.CNAs=np.array(MD.CNAs)
        savedata('cna',i, MD.CNAs)
    if MD.BOOon:
        MD.BOOs=np.array(MD.BOOs)
        savedata('boo',i, MD.BOOs)
    if MD.BOOtimeson:
        MD.BOOstimestemp=np.array(MD.BOOstimestemp)
        savedata('bootimes',i, MD.BOOstimestemp)
    if MD.VIon:
        #各列长度不同，无法保存成array
        savedata('vi',i, MD.VoronoiIndex)
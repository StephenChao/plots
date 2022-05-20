#!usr/bin/env python
import os
import glob
import math
import datetime
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse   import OptionParser
from time       import gmtime, strftime
from array import array
from ROOT import gROOT, TPaveLabel, TPie, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack, TGraph, TGraphErrors,TChain,TArrow, TCanvas, TMatrixDSym, TMath, TText, TPad, TVectorD, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite

parser = OptionParser()
parser.add_option('--channel',    action="store",type="string",dest="channel"    ,default="had")
parser.add_option('--MODE',       action="store",type="string",dest="MODE"       ,default="MC" )
parser.add_option('--REGION',     action="store",type="string",dest="REGION"     ,default="PS" )
parser.add_option('--TAG',        action="store",type="string",dest="TAG"        ,default=""   )
parser.add_option('--SFs',        action="store",type="int"   ,dest="SFs"        ,default=0    )
parser.add_option('--piechart',   action="store",type="int"   ,dest="piechart"   ,default=0    )
parser.add_option('--tau',        action="store",type="float" ,dest="tau"        ,default=0.4  )
parser.add_option('--y',          action="store",type="string",dest="y"          ,default="16,17,18")
parser.add_option('--FBT',        action="store",type="int"   ,dest="FBT"        ,default=0    )
(options, args) = parser.parse_args()

def UnderOverFlow1D(h):
    Bins=h.GetNbinsX();
    h.SetBinContent( 1,  h.GetBinContent(1)+h.GetBinContent(0) );  h.SetBinError(   1,  math.sqrt( h.GetBinError(1)*h.GetBinError(1) + h.GetBinError(0)*h.GetBinError(0)) );
    h.SetBinContent( Bins,  h.GetBinContent(Bins)+h.GetBinContent(Bins+1) );  h.SetBinError(   Bins,  math.sqrt( h.GetBinError(Bins)*h.GetBinError(Bins) + h.GetBinError(Bins+1)*h.GetBinError(Bins+1)) );
    return h;

def Integerization(h):
    Bins=h.GetNbinsX();
    for i in range(1,Bins+1):
        if (h.GetBinContent(i)-int(h.GetBinContent(i)))>0.5 : value=int(h.GetBinContent(i))+1;
        else : value=int(h.GetBinContent(i));
        h.SetBinContent( i, value          );
        h.SetBinError(   i, math.sqrt(value) );
    return h;

def Get_MeanVec_CovMat(S,   B ):
    S_Vec_of_Means=[];S_Cov_Matrix=[];B_Vec_of_Means=[];B_Cov_Matrix=[];

    S_Vec_of_Means.append( S.GetMean() );    S_Cov_Matrix.append( S.GetStdDev()*S.GetStdDev() );
    B_Vec_of_Means.append( B.GetMean() );    B_Cov_Matrix.append( B.GetStdDev()*B.GetStdDev() );

    print " ==>  S Means:",S_Vec_of_Means,"StDevs^2=Var:", S_Cov_Matrix;
    print " ==>  B Means:",B_Vec_of_Means,"StDevs^2=Var:", B_Cov_Matrix;

    return [ S_Vec_of_Means, S_Cov_Matrix, S_Vec_of_Means, S_Cov_Matrix ]; # List of 4 objects, 2 vectors, 2 matrices



def FBT_(total, h):
    Bins = h.GetNbinsX(); I = total.Integral();
    for i in range(1, Bins + 1):
        valTotal = total.GetBinContent(i); ErrorTotal = total.GetBinError(i); 
        if valTotal==0: continue;
        value = h.GetBinContent(i); Error = h.GetBinError(i);
        h.SetBinContent( i, value*(I/Bins)/valTotal );  h.SetBinError(   i, Error*(I/Bins)/valTotal );
    return H;

def SOT_1D(h,sig):
    Bins=h.GetNbinsX(); Bins_n_Sigs={}; h_clone = h.Clone("h_clone");
    for i in range(1, Bins+1) : Bins_n_Sigs[i]=sig.GetBinContent(i);
    k=1;
    for i in sorted(Bins_n_Sigs, key=Bins_n_Sigs.get) :
        h.SetBinContent(k, h_clone.GetBinContent(i) ); k=k+1;

def OptimalCut(B,S,NORM_s1):
    Bins = B.GetNbinsX(); B0=B.Integral(); S0=S.Integral(); SigMax=0; InitialSig=S0/((B0+1)**0.5); #LeftCutBin=1; RightCutBin=Bins; BKGrjc=0; Sig_Eff=0; SigMax_Print=0; S_o_sqrtB=0;
    for RightEnd in range(Bins,0,-1):
        for LeftEnd in range(1,RightEnd+1):            #print LeftEnd,RightEnd;
            sig = S.Integral( LeftEnd , RightEnd )/((B.Integral( LeftEnd , RightEnd )+1)**0.5);
            if sig>SigMax: 
                SigMax=sig; LeftCutBin=LeftEnd; RightCutBin=RightEnd; 
                BKGrjc  = float(round(100*(1-B.Integral(LeftCutBin,RightCutBin)/(B0+0.000001)),2) );
                Sig_Eff = float(round(100*   S.Integral(LeftCutBin,RightCutBin)/(S0+0.000001) ,1) );
                SigMax_Print=float(round(100*(SigMax-InitialSig)/InitialSig,1));
                S_o_sqrtB=float( round(SigMax/NORM_s1,3) );
    ResultList = [ LeftCutBin , RightCutBin , S.GetBinLowEdge(LeftCutBin) , S.GetBinLowEdge(RightCutBin+1), str(BKGrjc), str(Sig_Eff), str(SigMax_Print),str(S_o_sqrtB) ];   #print ResultList;
    #L_balance = (S.GetBinContent(LeftCutBin-1)+S.GetBinContent(LeftCutBin))  /(B.GetBinContent(LeftCutBin-1)+(B.GetBinContent(LeftCutBin)))   - 0.5*S.Integral(LeftCutBin , RightCutBin)/B.Integral(LeftCutBin , RightCutBin) ;
    #R_balance = (S.GetBinContent(RightCutBin-1)+S.GetBinContent(RightCutBin))/(B.GetBinContent(RightCutBin-1)+(B.GetBinContent(RightCutBin))) - 0.5*S.Integral(LeftCutBin , RightCutBin)/B.Integral(LeftCutBin , RightCutBin) ;
    #print "bin L",LeftCutBin,L_balance,"   -----    bin R",RightCutBin,R_balance; 
    return ResultList;  
    
def OptimalCut2(B,S,NORM_s1):
    Bins = B.GetNbinsX(); B0=B.Integral(); S0=S.Integral(); SigMax=0; LeftCutBin=0; RightCutBin=Bins; MidCutBin=Bins; InitialSig=S0/((B0+1)**0.5);
    for RightEnd in range(Bins,0,-1):
        for LeftEnd in range(1,RightEnd+1):
            for j in range(LeftEnd,RightEnd):
                sig_2bins = S.Integral( LeftEnd , j )/((B.Integral( LeftEnd , j )+1)**0.5) + S.Integral( j+1 , RightEnd )/((B.Integral( j+1 , RightEnd )+1)**0.5);
                if sig_2bins>SigMax:
                    SigMax=sig_2bins; LeftCutBin=LeftEnd; MidCutBin=j; RightCutBin=RightEnd;
    BKGrjc  = int(round(100*(1-B.Integral(LeftCutBin,RightCutBin)/(B0+0.000001))));
    Sig_Eff = int(round(100*   S.Integral(LeftCutBin,RightCutBin)/(S0+0.000001)) );
    SigMax_Print=int(round(100*(SigMax-InitialSig)/InitialSig));
    ResultList = [ LeftCutBin , RightCutBin , S.GetBinLowEdge(LeftCutBin) , S.GetBinLowEdge(RightCutBin+1), str(BKGrjc), str(Sig_Eff), str(SigMax_Print), S.GetBinLowEdge(MidCutBin+1) ];
    return ResultList;


def RationUnc(h_data,h_TotalMC,h_Ratio,MaxY):
    for i in range(1,h_Ratio.GetNbinsX()+1,1):
        D  = h_data.GetBinContent(i);    eD = h_data.GetBinError(i);
        if D==0: eD=0.92;
        B  = h_TotalMC.GetBinContent(i); eB = h_TotalMC.GetBinError(i);
        if B<0.1 and eB>=B : eB=0.92; Err= 0.;
        if B!=0.        :Err=TMath.Sqrt( (eD*eD)/(B*B)  +(D*D*eB*eB)/(B*B*B*B)     ); h_Ratio.SetBinContent(i, D/B   );  h_Ratio.SetBinError(i, Err); #print i,")",h_Ratio.GetNbinsX()+1,")   data:",D," pm ",eD,"     Bkg:",B," pm ",eB,"   R:",D/B," pm ", Err
        if B==0.        :Err=TMath.Sqrt( (eD*eD)/(eB*eB)+(D*D*eB*eB)/(eB*eB*eB*eB) ); h_Ratio.SetBinContent(i, D/0.92);  h_Ratio.SetBinError(i, Err);
        if D==0 and B==0:                                                             h_Ratio.SetBinContent(i, -1);      h_Ratio.SetBinError(i, 0  );
        if h_Ratio.GetBinContent(i)>MaxY:h_Ratio.SetBinContent(i,MaxY); ### To visualise the points above axis... #h_Ratio.Fit("pol1");
    return h_Ratio;



class ANALYSIS:
    def __init__(self, channel , fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1", input_workspace=None):
        self.setTDRStyle();
        self.channel    =channel;
        self.color_palet={'data':1, 'QCD':2,  'Rest':62,  'VV':62, 'STop':8, 'TTbar':80, 'ZJets':6, 'WJets':90, 'Signal':1, 'Uncertainty':1, }
        self.Signal_Scale1 =1;
        self.Signal_Scale2 =1;
        self.exponent = None;
        #self.Get_MeanVec_CovMat=[];
        # add by Qilong start
        self.KeepColumn = {"Pt_tag3","Pt_tag1", "MET_et", "Nj8", "weight"}
        self.Plot2DFaster = False
        self.Optimal = True
        # add by Qilong start

        # ========= VARIABLES =============
        self.Intime_Cut_Variable = ["DPhi_mgR_mg","DPhi_mg_MET","DPhi_mgR_MET","MET_o_PT_R","MR_v2"];
        # ========= Intime =============

    #================ SETTINGS FOR Canvas/pads/histos and more ==================
    def setTDRStyle(self):
        self.tdrStyle =TStyle("tdrStyle","Style for P-TDR");  self.tdrStyle.SetCanvasBorderMode(0);        self.tdrStyle.SetCanvasColor(kWhite);        self.tdrStyle.SetCanvasDefH(700);        self.tdrStyle.SetCanvasDefW(700);        self.tdrStyle.SetCanvasDefX(0);          self.tdrStyle.SetCanvasDefY(0);
        self.tdrStyle.SetPadBorderMode(0);             self.tdrStyle.SetPadColor(kWhite);        self.tdrStyle.SetPadGridX(False);        self.tdrStyle.SetPadGridY(False);        self.tdrStyle.SetGridColor(0);        self.tdrStyle.SetGridStyle(3);        self.tdrStyle.SetGridWidth(1);      
        self.tdrStyle.SetFrameBorderMode(0);        self.tdrStyle.SetFrameBorderSize(1);        self.tdrStyle.SetFrameFillColor(0);        self.tdrStyle.SetFrameFillStyle(0);        self.tdrStyle.SetFrameLineColor(1);        self.tdrStyle.SetFrameLineStyle(1);        self.tdrStyle.SetFrameLineWidth(1);
        self.tdrStyle.SetHistLineColor(1);        self.tdrStyle.SetHistLineStyle(0);        self.tdrStyle.SetHistLineWidth(1);        self.tdrStyle.SetEndErrorSize(2);              self.tdrStyle.SetMarkerStyle(20);      self.tdrStyle.SetErrorX(0.);
        self.tdrStyle.SetOptFit(1);        self.tdrStyle.SetFitFormat("5.4g");        self.tdrStyle.SetFuncColor(2);        self.tdrStyle.SetFuncStyle(1);        self.tdrStyle.SetFuncWidth(1);      self.tdrStyle.SetOptDate(0);      
        self.tdrStyle.SetOptFile(0); self.tdrStyle.SetOptStat(0); self.tdrStyle.SetStatColor(kWhite); self.tdrStyle.SetStatFont(42); self.tdrStyle.SetStatFontSize(0.025); self.tdrStyle.SetStatTextColor(1); self.tdrStyle.SetStatFormat("6.4g"); self.tdrStyle.SetStatBorderSize(1); self.tdrStyle.SetStatH(0.1); self.tdrStyle.SetStatW(0.15);
        self.tdrStyle.SetPadTopMargin(0.05);        self.tdrStyle.SetPadBottomMargin(0.13);        self.tdrStyle.SetPadLeftMargin(0.18);        self.tdrStyle.SetPadRightMargin(0.06);      
        self.tdrStyle.SetOptTitle(0);        self.tdrStyle.SetTitleFont(42);        self.tdrStyle.SetTitleColor(1);        self.tdrStyle.SetTitleTextColor(1);        self.tdrStyle.SetTitleFillColor(10);        self.tdrStyle.SetTitleFontSize(0.05);
        self.tdrStyle.SetTitleColor(1, "XYZ");        self.tdrStyle.SetTitleFont(42, "XYZ");        self.tdrStyle.SetTitleSize(0.06, "XYZ");  
        self.tdrStyle.SetTitleXOffset(0.8);        self.tdrStyle.SetTitleYOffset(0.8);      
        self.tdrStyle.SetLabelColor(1, "XYZ");        self.tdrStyle.SetLabelFont(42, "XYZ");        self.tdrStyle.SetLabelOffset(0.007, "XYZ");        self.tdrStyle.SetLabelSize(0.04, "XYZ");
        self.tdrStyle.SetAxisColor(1, "XYZ");        self.tdrStyle.SetStripDecimals(kTRUE);        self.tdrStyle.SetTickLength(0.03, "XYZ");        self.tdrStyle.SetNdivisions(510, "XYZ");        self.tdrStyle.SetPadTickX(1);       self.tdrStyle.SetPadTickY(1);      
        self.tdrStyle.SetOptLogx(0); self.tdrStyle.SetOptLogy(0); self.tdrStyle.SetOptLogz(0);
        self.tdrStyle.SetPaperSize(20.,20.); self.tdrStyle.cd();


    def DefineSelection_0lep(self):#==========[ 0lep REGIONS & SELECTION ]===========================================
        REGION=options.REGION; MODE=options.MODE; year=options.y;
        i=0.6;j=0.1;k=1.0;  
        deep_R ="("+str(i)+"*dnnDecorrw_a + "+str(j)+"*dnnDecorrh4q_a + "+str(k)+"*dnnDecorrtop_a) / (dnnDecorrqcd_a + "+str(i)+"*dnnDecorrw_a + "+str(j)+"*dnnDecorrh4q_a + "+str(k)+"*dnnDecorrtop_a)";
        i=0.03; j=0.03; k=0.9;  
        deep_g =" ("+str(i)+"*dnnDecorr_probQCDothers_c + "+str(j)+"*(dnnDecorr_probQCDb_c+dnnDecorr_probQCDc_c) + "+str(k)+"*(dnnDecorr_probQCDbb_c+dnnDecorr_probQCDcc_c)) / (dnnDecorr_probTbqq_c + ("+str(i)+"*dnnDecorr_probQCDothers_c + "+str(j)+"*(dnnDecorr_probQCDb_c+dnnDecorr_probQCDc_c) + "+str(k)+"*(dnnDecorr_probQCDbb_c+dnnDecorr_probQCDcc_c)) ) " 

        t41max = "tau41_max"; t41mid = "tau41_mid"; t41min = "tau41_min"; t31max = "tau31_max"; t31mid = "tau31_mid"; t31min = "tau31_min"; t21max = "tau21_max"; t21mid = "tau21_mid"; t21min = "tau21_min"; t42max = "tau42_max"; t42mid = "tau42_mid"; t42min = "tau42_min";
        Nb="nbtag_deep"; Nj4="num_ak4jetsex"; NbT="nbtag_deep_tight"; NbL="nbtag_deep_loose";       
        
        #deep_WH_max="(jetAK8puppi_dnnDecorrh4q_max+jetAK8puppi_dnnDecorrw_max)/(jetAK8puppi_dnnDecorrh4q_max+jetAK8puppi_dnnDecorrw_max+jetAK8puppi_dnnDecorrqcd_max)";        #deep_W_max ="jetAK8puppi_dnnDecorrW_max";   deep_W_mid ="jetAK8puppi_dnnDecorrW_mid";   deep_W_min ="jetAK8puppi_dnnDecorrW_min";        #deep_t_max ="jetAK8puppi_dnnDecorrTop_max"; deep_t_mid ="jetAK8puppi_dnnDecorrTop_mid"; deep_t_min ="jetAK8puppi_dnnDecorrTop_min";
        #deep_WH_a    ="(dnnDecorrh4q_a+dnnDecorrw_a)/(dnnDecorrqcd_a+dnnDecorrh4q_a+dnnDecorrw_a)";
        #deep_WHt_a   ="(dnnDecorrh4q_a+dnnDecorrw_a+dnnDecorrtop_a)/(dnnDecorrqcd_a+dnnDecorrh4q_a+dnnDecorrw_a+dnnDecorrtop_a)";
        #deep_Wt_a    ="(dnnDecorrw_a+dnnDecorrtop_a)/(dnnDecorrqcd_a+dnnDecorrw_a+dnnDecorrtop_a)";
        #deep_Ht_a    ="(dnnDecorrh4q_a+dnnDecorrtop_a)/(dnnDecorrqcd_a+dnnDecorrh4q_a+dnnDecorrtop_a)";
        #deep_MC_WH_a ="(dnnh4q_a+dnnw_a)/(dnnqcd_a+dnnh4q_a+dnnw_a)";
        #deep_MC_WHt_a="(dnnh4q_a+dnnw_a+dnntop_a)/(dnnqcd_a+dnnh4q_a+dnnw_a+dnntop_a)";
        #deep_Wt_a    ="(dnnw_a+dnntop_a)/(dnnqcd_a+dnnw_a+dnntop_a)";
        #deep_MC_Ht_a ="(dnnh4q_a+dnntop_a)/(dnnqcd_a+dnnh4q_a+dnntop_a)";
        #deep_WH_c    ="(dnnDecorrh4q_c+dnnDecorrw_c)/(dnnDecorrqcd_c+dnnDecorrh4q_c+dnnDecorrw_c)";
        #deep_WHt_c   ="(dnnDecorrh4q_c+dnnDecorrw_c+dnnDecorrtop_c)/(dnnDecorrqcd_c+dnnDecorrh4q_c+dnnDecorrw_c+dnnDecorrtop_c)";
        #deep_Wt_c    ="(dnnDecorrw_c+dnnDecorrtop_c)/(dnnDecorrqcd_c+dnnDecorrw_c+dnnDecorrtop_c)";
        #deep_Ht_c    ="(dnnDecorrh4q_c+dnnDecorrtop_c)/(dnnDecorrqcd_c+dnnDecorrh4q_c+dnnDecorrtop_c)";
        #deep_MC_WH_c ="(dnnh4q_c+dnnw_c)/(dnnqcd_c+dnnh4q_c+dnnw_c)";
        #deep_MC_WHt_c="(dnnh4q_c+dnnw_c+dnntop_c)/(dnnqcd_c+dnnh4q_c+dnnw_c+dnntop_c)";
        #deep_MC_Wt_c ="(dnnw_c+dnntop_c)/(dnnqcd_c+dnnw_c+dnntop_c)";
        #deep_MC_Ht_c ="(dnnh4q_c+dnntop_c)/(dnnqcd_c+dnnh4q_c+dnntop_c)";

        #-------------------------- PRE/SELECTION -----------------------------------------------------
        PS = "HT>1100 && ( (Nj8==2 && MJJ>900 ) || (Nj8==3 && MJJJ>900 ) || (Nj8==4 && MJJJ>900 ) ) && PTj>550 && PTj_2>200 && Mj_max>50 && nb_t_deep_in==0 && flag_Int==1 ";
        PS2="("+PS+") && Nj8==2 && Mj_a>50 && fabs(Etaj_a-Etaj_c)<2";
        PS3="("+PS+") && Nj8==3 && Mj_a>50 && Mj_b>50";
 
        Low_MET ="&& (                  MET_et/PTj_a<0.3)";
        Med_MET ="&& (0.3<MET_et/PTj_a&&MET_et/PTj_a<0.6)";
        High_MET="&& (0.6<MET_et/PTj_a                  )";

        SR1a=PS2+Low_MET+ "&&    BDT_R_MET_0_0p4  > .9 && BDT_g>.5 ";
        SR1b=PS2+Low_MET+ "&& ( (BDT_R_MET_0_0p4  > .9 && BDT_g<.5) || (0.75<BDT_R_MET_0_0p4   && BDT_R_MET_0_0p4  <0.9  && BDT_g>0.5) || (BDT_R_MET_0_0p4  <0.75 && BDT_g>0.8) ) ";
        SR2a=PS2+Med_MET+ "&&    BDT_R_MET_0p4_Inf>.85 && BDT_g>.5 ";
        SR2b=PS2+Med_MET+ "&& ( (BDT_R_MET_0p4_Inf>.85 && BDT_g<.5) || (0.6 <BDT_R_MET_0p4_Inf && BDT_R_MET_0p4_Inf<0.85 && BDT_g>0.5) || (BDT_R_MET_0p4_Inf<0.6  && BDT_g>0.8) ) ";
        SR3a=PS2+High_MET+"&&    BDT_R_MET_0p4_Inf>.85 && BDT_g>.5 ";
        SR3b=PS2+High_MET+"&& ( (BDT_R_MET_0p4_Inf>.85 && BDT_g<.5) || (0.6 <BDT_R_MET_0p4_Inf && BDT_R_MET_0p4_Inf<0.85 && BDT_g>0.5) || (BDT_R_MET_0p4_Inf<0.6  && BDT_g>0.8) ) ";
        SR4a=PS3+         "&&    BDT_g > 0.3 && PartNet_MD_W_a > 0.95 && PartNet_MD_W_b > 0.9  &&  60<Mj_a&&Mj_a<100  &&  60<Mj_b&&Mj_b<100 ";
        SR4b=PS3+         "&&    BDT_g > 0.3 && ( (PartNet_MD_W_a>.95&&PartNet_MD_W_b<.9) || (0.85<PartNet_MD_W_a&&PartNet_MD_W_a<.95 && 0.8<PartNet_MD_W_b) )  &&  60<Mj_a&&Mj_a<100  &&  60<Mj_b&&Mj_b<100 ";

        SR1ab="("+SR1a+") || ("+SR1b+")";
        SR2ab="("+SR2a+") || ("+SR2b+")";
        SR3ab="("+SR3a+") || ("+SR3b+")";
        SR4ab="("+SR4a+") || ("+SR4b+")";

        CR1 =PS2+Low_MET+ "&& ( (BDT_R_MET_0_0p4  < .9  && BDT_g<.5) || (BDT_R_MET_0_0p4  <0.75  && BDT_g<0.8) )";
        CR2 =PS2+Med_MET+ "&& ( (BDT_R_MET_0_0p4  < .85 && BDT_g<.5) || (BDT_R_MET_0_0p4  <0.6   && BDT_g<0.8) )";
        CR3 =PS2+High_MET+"&& ( (BDT_R_MET_0_0p4  < .85 && BDT_g<.5) || (BDT_R_MET_0_0p4  <0.6   && BDT_g<0.8) )";
        CR4 =PS3+         "&& (  BDT_g > 0.3 && !( (PartNet_MD_W_a>.95&&PartNet_MD_W_b<.9) || (0.85<PartNet_MD_W_a&&PartNet_MD_W_a<.95 && 0.8<PartNet_MD_W_b) )  &&  60<Mj_a&&Mj_a<100  &&  60<Mj_b&&Mj_b<100)";

        if MODE in ["MC","MCvsDATA","DECO"]: self.Make_Controlplots_for_0lep( eval(REGION),"","" );
        if MODE == "COMP"                  : self.Make_Controlplots_for_0lep( eval("C"+options.REGION[1:3]) ,eval("S"+options.REGION[1:]), "" );

        if MODE == "2DPlot":
            #for i in DNN_Scores:
            #    for j in DNN_Scores:
            #        i=i+"_a"; j=j+"_a";
            #        self.construct_2Dplot(eval(REGION), ("pow(dnnDecorr_probTbqq_a,0.169)"   ,"1-pow(dnnDecorr_probQCDothers_a,0.29)"), "Tbqq^0.169" , "1-QCDothers^0.29", 20,   0, 1, 20, 0, 1, Intime=True);
            #     if DNN_Names[j] in ["QCDothers","QCDc","QCDcc","QCDb","QCDbb","Tbc","Tbq"]:
            #            self.construct_plot(Nj, "1-pow("+i+","+str(Exponents[j])+")" ,selection,"","0_"+str(j+1), 40, 0, 0.8 ,"1-"+DNN_Names[j]+"^"+str(round(Exponents[j],3))+"    MD raw score" ,"Events",logy,CR,Intime=True);
            #        else:
            #            self.construct_plot(Nj, "pow("+i+","+str(Exponents[j])+")" ,selection,"","0_"+str(j+1), 40, 0, 0.8 , DNN_Names[j]+"^"+str(round(Exponents[j],3))+"    MD raw score" ,"Events",logy,CR,Intime=True);
            self.construct_2Dplot(eval(REGION), ("pow(dnnDecorr_probTbqq_a,0.169)"   ,"1-pow(dnnDecorr_probQCDothers_a,0.29)"), "Tbqq^0.169" , "1-QCDothers^0.29", 40,   0, 1, 40, 0, 1, Intime=True);
            self.construct_2Dplot(eval(REGION), ("pow(dnnDecorr_probTbqq_a,0.169)"   ,"1-pow(dnnDecorr_probQCDc_a,0.136)"), "Tbqq^0.169" , "1-QCDothers^0.136", 40,   0, 1, 40, 0, 1, Intime=True);

            #self.construct_2Dplot(eval(REGION), ("DPhi_j12"   ,"MET_et"), "DPhi_j12" , "MET_et", 30,   0, 3.2, 30, 0, 300 );
            #self.construct_2Dplot(eval(REGION), ("DR_j12"     ,"MET_et"), "DR_j12"   , "MET_et", 30,   0, 5  , 30, 0, 300 );
            #self.construct_2Dplot(eval(REGION), ("MET_o_PT_R" , "p20PT_g_o_PT_R" ), "MET / PT_R"         , "0.2 PT(g) / PT(R)"  , 20, .0,   .5, 20, 0,   1, Intime=True);
            #self.construct_2Dplot(eval(REGION), ("DPhi_mgR_R" , "DPhi_R_MET"   ), "|#Delta#Phi(R,-R-g)|" , "|#Delta#Phi(R,MET)|", 16,  0,  3.2, 15, 0, 3.2, Intime=True);
            #self.construct_2Dplot(eval(REGION), ("DPhi_mg_R"  , "DPhi_R_MET"   ), "|#Delta#Phi(R,-g)|"   , "|#Delta#Phi(R,MET)|", 16,  0,  3.2, 15, 0, 3.2, Intime=True);
            #self.construct_2Dplot(eval(REGION), ("PT_asy"        ,"MET_et")                         , "PT_asy", "MET_et", 30,   0, 300, 30, 0, 300, Intime=True);
            #self.construct_2Dplot(eval(REGION), ("PTj_c-PTj_a"   ,"MET_et")      , "PTj_c-PTj_a"    , "MET_et", 30,-300, 300, 30, 0, 300 );
            #self.construct_2Dplot(eval(REGION), ("PT_imb"        ,"MET_et")      , "PT_imb"         , "MET_et", 30,   0, 300, 30, 0, 300, Intime=True);
            #self.construct_2Dplot(eval(REGION), ("PT_asy_o_PTj_a","MET_et/PTj_a"), "PT_asy / PTj_a" , "MET_et / PTj_R", 30,   0, 6, 30, 0, 6, Intime=True);
            #self.construct_2Dplot(eval(REGION), ("PT_imb_o_PTj_a","MET_et/PTj_a"), "PT_imb / PTj_a" , "MET_et / PTj_R", 30,   0, 6, 30, 0, 6, Intime=True);
            # to plot faster (not recommended because this need more modification for the code for different variables)


    def Make_Controlplots_for_0lep(self,selection,selection2,tag,CR=0):
        REGION=options.REGION; Nj=234; MODE=options.MODE; tag=""; year=options.y; logy=0; IT=False; SEL=selection2;
        t41max = "tau41_max"; t41mid = "tau41_mid"; t41min = "tau41_min"; t31max = "tau31_max"; t31mid = "tau31_mid"; t31min = "tau31_min"; t21max = "tau21_max"; t21mid = "tau21_mid"; t21min = "tau21_min"; t42max = "tau42_max"; t42mid = "tau42_mid"; t42min = "tau42_min";
        Nb="nbtag_deep"; Nj4="num_ak4jetsex"; NbT="nbtag_deep_tight"; NbL="nbtag_deep_loose";          #deep_WH_max="(jetAK8puppi_dnnDecorrh4q_max+jetAK8puppi_dnnDecorrw_max)/(jetAK8puppi_dnnDecorrh4q_max+jetAK8puppi_dnnDecorrw_max+jetAK8puppi_dnnDecorrqcd_max)";        #deep_H_max ="jetAK8puppi_dnnDecorrH4q_max";deep_H_mid="jetAK8puppi_dnnDecorrH4q_mid";deep_H_min="jetAK8puppi_dnnDecorrH4q_min";        #deep_W_max ="jetAK8puppi_dnnDecorrW_max";  deep_W_mid="jetAK8puppi_dnnDecorrW_mid";  deep_W_min="jetAK8puppi_dnnDecorrW_min";        #deep_t_max ="jetAK8puppi_dnnDecorrTop_max";deep_t_mid="jetAK8puppi_dnnDecorrTop_mid";deep_t_min="jetAK8puppi_dnnDecorrTopmin";
        i=0.6;j=0.1;k=1.0;  #deep_R ="("+str(i)+"*dnnDecorrw_a + "+str(j)+"*dnnDecorrh4q_a + "+str(k)+"*dnnDecorrtop_a) / (dnnDecorrqcd_a + "+str(i)+"*dnnDecorrw_a + "+str(j)+"*dnnDecorrh4q_a + "+str(k)+"*dnnDecorrtop_a)";
        i=1; j=1; k=1;      deep_g ="("+str(i)+"*dnnDecorr_probQCDothers_c + "+str(j)+"*(dnnDecorr_probQCDb_c+dnnDecorr_probQCDc_c) + "+str(k)+"*(dnnDecorr_probQCDbb_c+dnnDecorr_probQCDcc_c)) / (dnnDecorr_probTbqq_c + ("+str(i)+"*dnnDecorr_probQCDothers_c + "+str(j)+"*(dnnDecorr_probQCDb_c+dnnDecorr_probQCDc_c) + "+str(k)+"*(dnnDecorr_probQCDbb_c+dnnDecorr_probQCDcc_c)) ) " 
        Intime=False;
        if MODE in ["MC", "MCvsDATA", "DECO", ]:
            if REGION[:3] in ["PS2","CR1","CR2","CR3","SR1","SR2","SR3"]:
                #self.construct_plot(Nj,"MJJ"       ,selection, "", tag,  26,1000,3600 ,"M(jj)  (GeV)", "Events" ,logy,CR);
                self.construct_plot(Nj,"MJJ_v1"    ,selection, "", tag  ,26,1000,3600 ,"M(jj)*  (GeV)" , "Events" ,logy,CR,Intime=True);
                #self.construct_plot(Nj,"Mj_a"      ,selection, "", tag  ,14,  50, 400 ,"M(ja)  (GeV)"  , "Events" ,logy,CR);
                #self.construct_plot(Nj,"MR_v1"     ,selection, "", tag  ,14,  50, 400 ,"M(R)*  (GeV)"  , "Events" ,logy,CR,Intime=True);
                #self.construct_plot(Nj,"Mj_c"      ,selection, "", tag  ,20,   0, 200 ,"M(jc)  (GeV)"   , "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_a"     ,selection, "", tag  ,24,  400,1600,"PT(ja)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_c"     ,selection, "", tag  ,24,  400,1600,"PT(jc)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj"       ,selection, "", tag  ,24,  400,1600,"PT(j1)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_3"     ,selection, "", tag  ,24,  400,1600,"PT(j2)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"fabs(Etaj_a-Etaj_c)",selection,selection2,tag,25,0,2.5,"|#Delta#etajj|", "Events" ,logy,CR);
                #self.construct_plot(Nj,"DPhi_g_R",selection,"",tag,32,0,3.2,"|#Delta#phijj|"          , "Events" ,logy,CR,Intime=True);
                #self.construct_plot(Nj,"X1"         ,selection,"",tag , 40 ,0.96, 1 , "X1"            , "Events" , logy,CR,Intime=True);
                #self.construct_plot(Nj,"BDT_R_MET_0_0p4",selection,"",tag, 20, 0, 1, "BDT for R->WW"  ,"Events", logy, CR);
            if REGION[:3] in ["PS3","CR4","SR4"]:
                #self.construct_plot(Nj,"MJJJ"      ,selection, "", tag ,28, 1000,3800, "M(jjj)  (GeV)", "Events" ,logy,CR);
                self.construct_plot(Nj,"MJJJ_v1"   ,selection, "", tag ,28, 1000,3800, "M(jjj)*  = M(jjj) -mja -mjb +160.8  (GeV)","Events" ,logy,CR,Intime=True);
                self.construct_plot(Nj,"MWW_v1"    ,selection, "", tag ,27,  300,3000, "M(WW)*  = M(jj) -mja -mjb +160.8  (GeV)","Events" ,logy,CR,Intime=True);
                #self.construct_plot(Nj,"Mj_a"      ,selection, "", tag ,14,   50, 400, "M(ja)  (GeV)" , "Events" ,logy,CR);
                #self.construct_plot(Nj,"Mj_b"      ,selection, "", tag ,14,   50, 400, "M(ja)  (GeV)" , "Events" ,logy,CR);
                #self.construct_plot(Nj,"Mj_c"      ,selection, "", tag ,20,    0, 200, "M(jc)  (GeV)" , "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_a"     ,selection, "", tag ,24,  400,1600, "PT(ja)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_c"     ,selection, "", tag ,24,  400,1600, "PT(jc)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj"       ,selection, "", tag ,24,  400,1600, "PT(j1)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_3"     ,selection, "", tag ,24,  400,1600, "PT(j2)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"fabs(Etaj_a-Etaj_c)",selection,"",tag,25,0,2.5,"|#Delta#etajj|","Events",logy,CR);
        if MODE in [ "COMP" ]:
            if REGION[:3] in ["PS2","CR1","CR2","CR3","SR1","SR2","SR3"]:
                #self.construct_plot(Nj,"MJJ"       ,selection, selection2, tag,  26,1000,3600 ,"M(jj)  (GeV)", "Events" ,logy,CR);
                self.construct_plot(Nj,"MJJ_v1"    ,selection, selection2, tag  ,26,1000,3600 ,"M(jj)*  (GeV)" , "Events" ,logy,CR,Intime=True);
                #self.construct_plot(Nj,"Mj_a"      ,selection, selection2, tag  ,14,  50, 400 ,"M(ja)  (GeV)"  , "Events" ,logy,CR);
                #self.construct_plot(Nj,"MR_v1"     ,selection, selection2, tag  ,14,  50, 400 ,"M(R)*  (GeV)"  , "Events" ,logy,CR,Intime=True);
                #self.construct_plot(Nj,"Mj_c"      ,selection, selection2, tag  ,20,   0, 200 ,"M(jc)  (GeV)"   , "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_a"     ,selection, selection2, tag  ,24,  400,1600,"PT(ja)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_c"     ,selection, selection2, tag  ,24,  400,1600,"PT(jc)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj"       ,selection, selection2, tag  ,24,  400,1600,"PT(j1)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_3"     ,selection, selection2, tag  ,24,  400,1600,"PT(j2)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"fabs(Etaj_a-Etaj_c)",selection,selection2,tag,25,0,2.5,"|#Delta#etajj|", "Events" ,logy,CR);
                #self.construct_plot(Nj,"DPhi_g_R",selection,selection2,tag,32,0,3.2,"|#Delta#phijj|"          , "Events" ,logy,CR,Intime=True);
                #self.construct_plot(Nj,"X1"         ,selection,selection2,tag , 40 ,0.96, 1 , "X1"            , "Events" , logy,CR,Intime=True);
                #self.construct_plot(Nj,"BDT_R_MET_0_0p4",selection,selection2,tag, 20, 0, 1, "BDT for R->WW"  ,"Events", logy, CR);
            if REGION[:3] in ["PS3","CR4","SR4"]:
                #self.construct_plot(Nj,"MJJJ"      ,selection, selection2, tag ,28, 1000,3800, "M(jjj)  (GeV)", "Events" ,logy,CR);
                self.construct_plot(Nj,"MJJJ_v1"   ,selection, selection2, tag ,28, 1000,3800, "M(jjj)*  = M(jjj) -mja -mjb +160.8  (GeV)","Events" ,logy,CR,Intime=True);
                self.construct_plot(Nj,"MWW_v1"    ,selection, selection2, tag ,27,  300,3000, "M(WW)*  = M(jj) -mja -mjb +160.8  (GeV)","Events" ,logy,CR,Intime=True);
                #self.construct_plot(Nj,"Mj_a"      ,selection, selection2, tag ,14,   50, 400, "M(ja)  (GeV)" , "Events" ,logy,CR);
                #self.construct_plot(Nj,"Mj_b"      ,selection, selection2, tag ,14,   50, 400, "M(ja)  (GeV)" , "Events" ,logy,CR);
                #self.construct_plot(Nj,"Mj_c"      ,selection, selection2, tag ,20,    0, 200, "M(jc)  (GeV)" , "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_a"     ,selection, selection2, tag ,24,  400,1600, "PT(ja)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_c"     ,selection, selection2, tag ,24,  400,1600, "PT(jc)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj"       ,selection, selection2, tag ,24,  400,1600, "PT(j1)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"PTj_3"     ,selection, selection2, tag ,24,  400,1600, "PT(j2)  (GeV)", "Events" ,logy,CR);
                #self.construct_plot(Nj,"fabs(Etaj_a-Etaj_c)",selection,selection2,tag,25,0,2.5,"|#Delta#etajj|","Events",logy,CR);

    def construct_2Dplot(self, cut, variable,Xtitle,Ytitle,Xnbin,Xmin,Xmax,Ynbin,Ymin,Ymax, KeepColumn = None, Intime = False):
        year=options.y;
        #2DP
        PT_imb = '''bool FLIP = (Mj_a<50 && (dnnDecorrW_a-dnnDecorrW_c)<0.4);  TLorentzVector g, R, MET, P_imb, P_asy;  if(FLIP){ R.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); g.SetPtEtaPhiM( PTj_a, Etaj_a, Phij_a, Mj_a ); } else{ R.SetPtEtaPhiM( PTj_a, Etaj_a, Phij_a, Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); } MET.SetPtEtaPhiM( MET_et, R.Eta(), MET_phi, 0 ); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
        return (P_imb).Pt(); '''
        PT_asy = '''bool FLIP = (Mj_a<50 && (dnnDecorrW_a-dnnDecorrW_c)<0.4);  TLorentzVector g, R, MET, P_imb, P_asy;  if(FLIP){ R.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); g.SetPtEtaPhiM( PTj_a, Etaj_a, Phij_a, Mj_a ); } else{ R.SetPtEtaPhiM( PTj_a, Etaj_a, Phij_a, Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); } MET.SetPtEtaPhiM( MET_et, R.Eta(), MET_phi, 0 ); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),R.Eta(),(-g-R).Phi(),0);
        return (-g-R).Pt(); '''
        PT_imb_o_PTj_a = '''bool FLIP = (Mj_a<50 && (dnnDecorrW_a-dnnDecorrW_c)<0.4);  TLorentzVector g, R, MET, P_imb, P_asy;  if(FLIP){ R.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); g.SetPtEtaPhiM( PTj_a, Etaj_a, Phij_a, Mj_a ); } else{ R.SetPtEtaPhiM( PTj_a, Etaj_a, Phij_a, Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); } MET.SetPtEtaPhiM( MET_et, R.Eta(), MET_phi, 0 ); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
        return P_imb.Pt() / R.Pt() ; '''
        PT_asy_o_PTj_a = '''bool FLIP = (Mj_a<50 && (dnnDecorrW_a-dnnDecorrW_c)<0.4);  TLorentzVector g, R, MET, P_imb, P_asy;  if(FLIP){ R.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); g.SetPtEtaPhiM( PTj_a, Etaj_a, Phij_a, Mj_a ); } else{ R.SetPtEtaPhiM( PTj_a, Etaj_a, Phij_a, Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); } MET.SetPtEtaPhiM( MET_et, R.Eta(), MET_phi, 0 ); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
        return P_asy.Pt() / R.Pt() ; '''
        MET_o_PT_R = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
        return  MET.Pt()/R.Pt() ; '''

        CircularSignificance = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
        return  0.5*sqrt( (MET.Pt()/R.Pt())*(MET.Pt()/R.Pt()) + (0.20*g.Pt()/R.Pt())*(0.20*g.Pt()/R.Pt()) ) ; '''

        MET_o_PT_R = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
        return  MET.Pt()/R.Pt() ; '''

        Bin = ("Plot2D",";%s;%s"%(Xtitle,Ytitle),Xnbin,Xmin,Xmax,Ynbin,Ymin,Ymax)

        # first parameter of each element need to be increase
        ContourLevel = [(0.1,ROOT.kWhite),(0.333,ROOT.kYellow),(0.667,ROOT.kOrange),(0.9,ROOT.kRed)]
        # sort the ContourLevel in case this is not increase
        ContourLevel.sort(key = lambda a: a[0], reverse=False)        

        REGION = options.REGION

        #if options.y=="17"      :path="/eos/user/y/yusong/qilong/NTuple_Output/17/all/Tree_9_24/";  lumi= 41.5;
        #if options.y=="16"      :path="/eos/user/y/yusong/qilong/ROOTFILE/PlotTree/V1/";            lumi= 36.3;
        #if options.y=="16"      :path="/afs/cern.ch/work/q/qiguo/NTuple/V2/";                       lumi= 36.3;
        if options.y=="16"      :path="/eos/user/q/qiguo/ROOTFILES/GKK/0lepton/Tree_BDT/UL16/";     lumi= 36.3;   #15/4/22
        #if options.y=="17"      :path="/eos/user/y/yusong/qilong/NTuple_Output/17/all/Tree/";  lumi= 41.5;
        #  For HWW
        #if options.y=="18"      :path="/eos/user/y/yuzhe/HWW/frame/Custom/Tree/newTree/";  lumi=59.81;

        # prepare Input File
        BKGs = ROOT.std.vector('string')()
        Files = [
            "PKUTree_QCD_"+year+".root"
            "PKUTree_Rest_"+year+".root"
            "PKUTree_ST_"+year+".root"
            "PKUTree_TT_"+year+".root"
            "PKUTree_WJetsToQQ_"+year+".root"
            ]

        F_Signal1 = path+"Res1ToRes2GluToGluVV_M1-2000_R-200_"+year+".root";
        F_Signal2 = path+"Res1ToRes2GluToGluVV_M1-3000_R-1500_"+year+".root";
        Files     = [path+i for i in Files]

        for n in Files:
            BKGs.push_back(n)
            if KeepColumn:
                dfBKG = ROOT.RDataFrame("PKUTree", BKGs, KeepColumn)
            else:
                dfBKG = ROOT.RDataFrame("PKUTree", BKGs)
        
        dfBKG = ROOT.RDataFrame("PKUTree", BKGs      )
        df    = ROOT.RDataFrame("PKUTree", F_Signal1 )
        df    = ROOT.RDataFrame("PKUTree", F_Signal2 )    
        
        for n in Files:
            BKGs.push_back(n)

        for Added_var in self.Intime_Cut_Variable:
            if Added_var in cut:
                if Added_var not in [str(i) for i in dfBKG.GetColumnNames()]:
                    dfBKG     = dfBKG.Define(Added_var  , eval(Added_var))
                    df        = df.Define(Added_var, eval(Added_var))


        # apply selection
        print "selection ==>", cut
        dfBKG = dfBKG.Filter(cut)
        df    = df.Filter(cut)

        Intime_plot_var = list(variable)
        if variable[0] not in [str(i) for i in dfBKG.GetColumnNames()]:
            Intime_plot_var[0] = "Var1"
            try:
                print eval(variable[0])
                dfBKG  = dfBKG.Define(    "Var1",eval(variable[0]))
                df     = df.Define(    "Var1",eval(variable[0]))
            except (NameError,SyntaxError):
                dfBKG  = dfBKG.Define(    "Var1","return (%s);"%(variable[0]))
                df     = df.Define(    "Var1","return (%s);"%(variable[0]))

        if variable[1] not in [str(i) for i in dfBKG.GetColumnNames()]:
            Intime_plot_var[1] = "Var2"
            try:
                print eval(variable[1])
                dfBKG  = dfBKG.Define("Var2",eval(variable[1]))
                df     = df.Define(    "Var2",eval(variable[1]))
            except (NameError,SyntaxError):
                dfBKG  = dfBKG.Define(    "Var2","return (%s);"%(variable[1]))
                df     = df.Define(    "Var2","return (%s);"%(variable[1]))

        canvas = ROOT.TCanvas("Canvas","Canvas", 700,700);

        h_Signal  =    df.Histo2D( Bin, Intime_plot_var[0], Intime_plot_var[1], "weight"); 
        h_TotalMC = dfBKG.Histo2D( Bin, Intime_plot_var[0], Intime_plot_var[1], "weight"); 
        #h_SoB     =    df.Histo2D( Bin, Intime_plot_var[0], Intime_plot_var[1], "weight"); h_SoB.Divide(h_TotalMC);

        #print "h_Signal, h_TotalMC =",h_Signal, h_TotalMC;

        h_Signal.SetStats(0); h_TotalMC.SetStats(0);
        h_Signal.Smooth();
            
        h_Signal.SetContour(1);
        maxZ = h_Signal.GetMaximum();
        h_List = [];

        for i,j in enumerate(ContourLevel):
            exec('hc_{i} = h_Signal.Clone("hc_{i}")'.format(i = i));
            exec('hc_{i}.SetContourLevel(0,maxZ*ContourLevel[i][0])'.format(i = i));
            exec('hc_{i}.SetLineColor(ContourLevel[i][1])'.format(i = i));
            exec('h_List.append(hc_{i})'.format(i = i));

        # Draw S and B superimposed
        h_TotalMC.Draw("COLZ");

        #  Draw S/B ratio:
        #h_SoB.Draw("COLZ");

        for i in h_List: i.Draw("same CONT3");

        variable1=variable[0];variable2=variable[1]
        for c in [".","/","(",")","[","]","*","+",">","<"," ","=",",",":","deep","dnn","Decorr","jetAK8puppi","ass_tag","t_tag","_tag","TMath","Cos","Sin"]:
            variable1=variable1.replace(c,"_")
            variable2=variable2.replace(c,"_")
        for c in ["__","___","____","_____","______","_"]:
            variable1=variable1.replace(c,"");
            variable2=variable2.replace(c,"");
        Name=variable1+"_"+REGION+"__"+variable2+".png"

        canvas.SaveAs(Name)
        print "\n --> generate %s &"%(Name);

    # add by Qilong end
    
    



    def construct_plot(self,Nj,variable,cut,cut1,tag,nbin,min,max,xtitle="",ytitle="",logy=1,CR=0, Intime = False):
        
        SFs=options.SFs; channel=options.channel; MODE=options.MODE;  REGION=options.REGION;  year=options.y;
        print "\n -->  MODE:",MODE," variable:",variable,"\n       { "+cut+" }\n";

        #----------------- paths to root files -------------------
        #if options.y=="16"      :path="/eos/cms/store/user/.........;  lumi=36.3;
        # if options.y=="17"      :path="/eos/user/y/yusong/qilong/NTuple_Output/17/all/Tree/mu";  lumi= 41.5;
        #if options.y=="17"      :path="/eos/user/y/yusong/qilong/NTuple_Output/17/all/Tree_9_24/mu";  lumi= 41.5;
        if options.y=="16"      :path="/eos/user/q/qiguo/ROOTFILES/GKK/0lepton/Tree_BDT/UL16/";            lumi= 36.3;
        #if options.y=="18"      :path="/eos/cms/store/user/.........;  lumi=59.7;
        #if options.y=="16,17,18":path="/eos/cms/store/user/.........;  lumi=138;

        #====== DEFINE CANVAS ==========================
        if MODE in ["MC","MCvsDATA","COMP"]:
            print " ---> Define canvas"; 
            canvas_controlplot = TCanvas(REGION+"_"+variable, REGION+"_"+variable, 700,700);
            fPads1 = TPad("pad1", "", 0.0, 0.29, 1.00, 1.00);
            fPads2 = TPad("pad2", "", 0.0, 0.00, 1.00, 0.29);
            fPads1.SetBottomMargin(0.007);fPads1.SetLeftMargin( 0.10);fPads1.SetRightMargin( 0.03);
            fPads2.SetLeftMargin(  0.10 );fPads2.SetRightMargin(0.03);fPads2.SetBottomMargin(0.25);
            fPads1.Draw(); fPads2.Draw(); fPads1.cd();
        if MODE in ["DECO"]:
            canvas_controlplot = TCanvas(REGION+"_"+variable, REGION+"_"+variable, 700,565);
            canvas_controlplot.SetLeftMargin(0.1); canvas_controlplot.SetRightMargin(0.03);

        #====================== DEFINE TREES AND HISTOS ======================================
        if MODE in ["MCvsDATA","COMP"]:
            print " ---> tree for data"; 
            t_data  = TChain("PKUTree"); t_data.Add(path+"PKUTree_data_"+year+".root");
        h_data    = TH1D("h_data"  ,"h_DATA"  +";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_data.Sumw2();  #c_data=h_data.Clone("c_data");
        t_Signal1 = TChain("PKUTree"); t_Signal1.Add(path+"PKUTree_gKK-2000_R-200_" +year+".root");  h_Signal1=TH1D("h_Signal1","h_Signal1"+";%s;%s"%(xtitle,ytitle),nbin,min,max);  h_Signal1.Sumw2();
        t_Signal2 = TChain("PKUTree"); t_Signal2.Add(path+"PKUTree_gKK-3000_R-1500_"+year+".root");  h_Signal2=TH1D("h_Signal2", "h_Signal2" + ";%s;%s" % (xtitle, ytitle), nbin, min, max); h_Signal2.Sumw2();
        t_QCD     = TChain("PKUTree"); t_QCD.Add(    path+"PKUTree_QCD_"+year+".root"            ); h_QCD  =TH1D("h_QCD"  ,"h_QCD"  +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_QCD.Sumw2();   
        t_WJets   = TChain("PKUTree"); t_WJets.Add(  path+"PKUTree_WJetsToQQ_"+year+".root"      ); h_WJets=TH1D("h_WJets","h_WJets"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_WJets.Sumw2(); 
        t_TTbar   = TChain("PKUTree"); t_TTbar.Add(  path+"PKUTree_TT_"+year+".root"             ); h_TTbar=TH1D("h_TTbar","h_TTbar"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_TTbar.Sumw2(); 
        t_STop    = TChain("PKUTree"); t_STop.Add(   path+"PKUTree_ST_"+year+".root"             ); h_STop =TH1D("h_STop" ,"h_STop" +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_STop.Sumw2();  
        t_Rest    = TChain("PKUTree"); t_Rest.Add(   path+"PKUTree_Rest_"+year+".root"           ); h_Rest =TH1D("h_Rest" ,"h_Rest" +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_Rest.Sumw2();  

        if Intime:
            # ROOT.EnableImplicitMT() # allow to use mutiple core
            MJJ_v1 = ''' TLorentzVector g,R,MET,P_imb;  R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c ); MET.SetPtEtaPhiM( MET_et,0,MET_phi,0 );  if( fabs(MET.DeltaPhi(-g)) < 1  &&  (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM( MET_et,R.Eta(),MET_phi,0 ); R=R+MET; }; 
                        P_imb.SetPtEtaPhiM( (-R-g-MET ).Pt(),0,(-R-g-MET).Phi() ,0 );
                        return ( g + R ).M(); '''
            MJJ_v2 = ''' TLorentzVector g,R,MET,P_imb;  R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c ); MET.SetPtEtaPhiM( MET_et,0,MET_phi,0 );  if( fabs(MET.DeltaPhi(-g)) < 1  &&  (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM( MET_et,R.Eta(),MET_phi,0 ); R=R+MET; }; 
                        P_imb.SetPtEtaPhiM( (-R-g-MET ).Pt(),0,(-R-g-MET).Phi() ,0 );
                        return ( g + R ).M() - g.M(); '''
            MR_v1  = ''' TLorentzVector g,R,MET,P_imb;  R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c ); MET.SetPtEtaPhiM( MET_et,R.Eta(),MET_phi,0 ); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),0,(-R-MET-g).Phi(),0 );
                        if( fabs(MET.DeltaPhi(-g)) < 1  &&  (MET.Pt()/R.Pt()) > 0.1 ){ R=R+MET; };
                        return ( R ).M(); '''

            MJJJ_v1 =''' TLorentzVector g,Wa,Wb,P_imb;  Wa.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); Wb.SetPtEtaPhiM( PTj_b,Etaj_b,Phij_b,Mj_b ); g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c );  return ( (g+Wa+Wb).M() -Wa.M()-Wb.M()+160.8 ); '''
            MJJJ_v2 =''' TLorentzVector g,Wa,Wb,P_imb;  Wa.SetPtEtaPhiM( PTj_a*(80.4/Mj_a),Etaj_a,Phij_a,80.4); Wb.SetPtEtaPhiM(PTj_b*(80.4/Mj_b),Etaj_b,Phij_b,80.4);  g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c );                 return (g+Wa+Wb).M(); '''
            MJJJ_v3 =''' TLorentzVector g,Wa,Wb,P_imb;  Wa.SetPtEtaPhiM( PTj_a*(80.4/Mj_a),Etaj_a,Phij_a,80.4); Wb.SetPtEtaPhiM(PTj_b*(80.4/Mj_b),Etaj_b,Phij_b,80.4);  g.SetPtEtaPhiM( PTj_c,Etaj_c,(-Wa-Wb).Phi(),Mj_c );         return (g+Wa+Wb).M(); '''
            MJJJ_v4 =''' TLorentzVector g,Wa,Wb,P_imb;  Wa.SetPtEtaPhiM( PTj_a*(80.4/Mj_a),Etaj_a,Phij_a,80.4); Wb.SetPtEtaPhiM(PTj_b*(80.4/Mj_b),Etaj_b,Phij_b,80.4);  g.SetPtEtaPhiM( (-Wa-Wb).Pt(),Etaj_c,(-Wa-Wb).Phi(),Mj_c ); return (g+Wa+Wb).M(); '''

            MWW_v1 =''' TLorentzVector g,Wa,Wb,P_imb;  Wa.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); Wb.SetPtEtaPhiM( PTj_b,Etaj_b,Phij_b,Mj_b ); g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c ); P_imb.SetPtEtaPhiM( (-g-Wa-Wb).Pt(),(-Wa-Wb).Eta(),(-g-Wa-Wb).Phi(),0 );  return ((Wa+Wb).M() -Wa.M()-Wb.M()+160.8 ); '''
            MWW_v2 =''' TLorentzVector g,Wa,Wb,P_imb;  Wa.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); Wb.SetPtEtaPhiM( PTj_b,Etaj_b,Phij_b,Mj_b ); g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c ); P_imb.SetPtEtaPhiM( (-g-Wa-Wb).Pt(),(-Wa-Wb).Eta(),(-g-Wa-Wb).Phi(),0 );  return  (Wa+Wb).M(); '''
            MWW_v3 =''' TLorentzVector g,Wa,Wb,P_imb;  Wa.SetPtEtaPhiM( PTj_a*(80.4/Mj_a),Etaj_a,Phij_a, 80.4 ); Wb.SetPtEtaPhiM( PTj_b*(80.4/Mj_b),Etaj_b,Phij_b,80.4 );                                                 g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c );  return  (Wa+Wb).M(); '''
            MWW_v4 =''' TLorentzVector g,Wa,Wb,P_imb;  Wa.SetPtEtaPhiM( PTj_a*(80.4/Mj_a),Etaj_a,Phij_a, 80.4 ); Wb.SetPtEtaPhiM( PTj_b,Etaj_b,Phij_b,Mj_b );                                                             g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c );  return  (Wa+Wb).M(); '''
            MWW_v5 =''' TLorentzVector g,Wa,Wb,P_imb;  Wa.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); Wb.SetPtEtaPhiM( PTj_b,Etaj_b,Phij_b,Mj_b ); g.SetPtEtaPhiM( PTj_c,Etaj_c,Phij_c,Mj_c ); P_imb.SetPtEtaPhiM( (-g-Wa-Wb).Pt(),(-Wa-Wb).Eta(),(-g-Wa-Wb).Phi(),0 );  return ((Wa+Wb).M() -Wa.M()+80.4 ); '''

            MTR   = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (R).Mt(); '''
            MTR2  = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);            return (R+MET).Mt(); '''

            DPhi_g_R   = ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0);
            return fabs((g).DeltaPhi(R)); '''
            DPhi_mg_R   = ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0);
            return fabs((-g).DeltaPhi(R)); '''
            DPhi_mg_MET = ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0);
            return fabs((-g).DeltaPhi(MET)); '''
            DPhi_mgR_MET = ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0);
            return fabs((-g-R).DeltaPhi(MET)); '''
            DPhi_mgR_Pim = ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0);
            return fabs((-g-R).DeltaPhi(P_imb)); '''
            DPhi_mgR_R= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0);
            return fabs((-g-R).DeltaPhi(R)); '''
            DPhi_R_MET= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0);
            return fabs((MET).DeltaPhi(R)); '''
            DPhi_R_Pim   = ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0);
            return fabs((R).DeltaPhi(P_imb)); '''
            DPhi_MET_Pim   = ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0);
            return fabs((MET).DeltaPhi(P_imb)); '''

            X1= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return    sqrt( 2*PTj_a*PTj_c* ( cosh(fabs(Etaj_a-Etaj_c))-cos(fabs(R.DeltaPhi(g))) ) ) / (R+g).M()    ; '''
            X2= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return    sqrt( 2*PTj_a*PTj_c* ( cosh(fabs(Etaj_a-Etaj_c))+1                        ) ) / (R+g).M()    ; '''

            X3= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return    sqrt(PTj_a*PTj_c) / min(PTj_a,PTj_c)   ; '''
            X4= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return    sqrt( 2*PTj_a*PTj_c* ( cosh(fabs(Etaj_a-Etaj_c))-cos(fabs(R.DeltaPhi(g))) ) )  /  sqrt( (R+g).M()*(R+g).M() -Mj_a*Mj_a - Mj_c*Mj_c )  ; '''

            cosh_cos= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  cosh(fabs(Etaj_a-Etaj_c))-cos(fabs(R.DeltaPhi(g)))  ; '''
            cosh= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  cosh(fabs(Etaj_a-Etaj_c)); '''
            Dphi_R_g= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  fabs(R.DeltaPhi(g))  ; '''
            PT_o_MJJ= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  sqrt( PTj_a*PTj_c ) / (g+R).M() ; '''
            PTmax_o_MJJ= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  max(PTj_a,PTj_c) / (g+R).M() ; '''
            PTmin_o_MJJ= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  min(PTj_a,PTj_c) / (g+R).M() ; '''
            X1min= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  sqrt( 2*min(PTj_a,PTj_c)*min(PTj_a,PTj_c)* ( cosh(fabs(Etaj_a-Etaj_c))-cos(fabs(R.DeltaPhi(g))) ) ) / (R+g).M()    ; '''
            X1max=''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  sqrt( 2*max(PTj_a,PTj_c)*max(PTj_a,PTj_c)* ( cosh(fabs(Etaj_a-Etaj_c))-cos(fabs(R.DeltaPhi(g))) ) ) / (R+g).M()    ; '''
            X2min= ''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  sqrt( 2*min(PTj_a,PTj_c)*min(PTj_a,PTj_c)* ( cosh(fabs(Etaj_a-Etaj_c))+1 ) ) / (R+g).M() ; '''
            X2max=''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return  sqrt( 2*max(PTj_a,PTj_c)*max(PTj_a,PTj_c)* ( cosh(fabs(Etaj_a-Etaj_c))+1 ) ) / (R+g).M() ; '''

            PT_asym=''' TLorentzVector g,R,MET,P_imb,P_asy;  R.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a);  g.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c);  MET.SetPtEtaPhiM(MET_et,0,MET_phi,0); if( fabs(MET.DeltaPhi(-g)) < 1 && (MET.Pt()/R.Pt()) > 0.1 ){ MET.SetPtEtaPhiM(MET_et,Etaj_a,MET_phi,0); };  P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),(-R-MET-g).Eta(),(-R-MET-g).Phi(),0 );  P_asy.SetPtEtaPhiM((-g-R).Pt(),(-R-g).Eta(),(-g-R).Phi(),0); 
            return PTj_3/PTj_1 ; '''

            DPhi_g_MET = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return fabs((g).DeltaPhi(MET)); '''
            DPhi_g_MET_v2 = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return fabs((-g).DeltaPhi(MET)); '''
            DR_g_R = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (g).DeltaR(R); '''
            DR_R_MET = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (R).DeltaR(MET); '''

            DR_R_P_imb = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return R.DeltaR(P_imb); '''
            DR_R_Pasy = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return R.DeltaR(P_asy); '''
            DR_mg_Pasy = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (-g).DeltaR(P_asy); '''
            DR_mg_R = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (-g).DeltaR(R); '''
            DPhi_PT_asy_R = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return fabs((R).DeltaPhi(P_asy)); '''
            PT_imb = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (P_imb).Pt(); '''
            PT_asy = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (P_asy).Pt(); '''
            PT_imb_o_PTj_a = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return P_imb.Pt() / R.Pt() ; '''
            PT_asy_o_PTj_a = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return P_asy.Pt() / R.Pt() ; '''
            PT_dif_gR = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (g.Pt()-R.Pt()); '''
            PT_g_R_MET = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return ( g.Pt() - R.Pt() + MET.Pt() ); '''
            Dot_R_P_imb = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (  R.Dot(P_imb)/(R.Mag()*P_imb.Mag()) ); '''
            Dot_R_MET = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return (  R.Dot(MET)/(R.Mag()*MET.Mag())  ); '''
            MET_o_PT_R = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return  MET.Pt()/R.Pt() ; '''
            p20PT_g_o_PT_R = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return  0.20*g.Pt()/R.Pt() ; '''
            CircularSignificance = '''bool FLIP = 0;  TLorentzVector g,R,MET,P_imb,P_asy; if(FLIP){ R.SetPtEtaPhiM(PTj_c,Etaj_c,Phij_c,Mj_c); g.SetPtEtaPhiM(PTj_a,Etaj_a,Phij_a,Mj_a); } else{ R.SetPtEtaPhiM( PTj_a,Etaj_a,Phij_a,Mj_a ); g.SetPtEtaPhiM( PTj_c, Etaj_c, Phij_c , Mj_c ); }; MET.SetPtEtaPhiM(MET_et,R.Eta(),MET_phi,0); P_imb.SetPtEtaPhiM( (-R-MET-g).Pt(),R.Eta(),(-R-MET-g).Phi(),0 ); P_asy.SetPtEtaPhiM((-g-R).Pt(),(-g-R).Eta(),(-g-R).Phi(),0);
            return  0.5*sqrt( (MET.Pt()/R.Pt())*(MET.Pt()/R.Pt()) + (0.20*g.Pt()/R.Pt())*(0.20*g.Pt()/R.Pt()) ) ; '''

            if MODE in ["MCvsDATA"] and REGION[:2]=="SR" : MODE="MC";  ## BLINDING

            if MODE in ["MCvsDATA","COMP"]:
                print " ---> tree for data";
                df_data = ROOT.RDataFrame("PKUTree", path+"PKUTree_data_"           +year+".root");
            df_QCD      = ROOT.RDataFrame("PKUTree", path+"PKUTree_QCD_"            +year+".root");
            df_WJets    = ROOT.RDataFrame("PKUTree", path+"PKUTree_WJetsToQQ_"      +year+".root");
            df_TTbar    = ROOT.RDataFrame("PKUTree", path+"PKUTree_TT_"             +year+".root");
            df_STop     = ROOT.RDataFrame("PKUTree", path+"PKUTree_ST_"             +year+".root");
            df_Rest     = ROOT.RDataFrame("PKUTree", path+"PKUTree_Rest_"           +year+".root");
            df_Signal1  = ROOT.RDataFrame("PKUTree", path+"PKUTree_gKK-2000_R-200_" +year+".root");
            df_Signal2  = ROOT.RDataFrame("PKUTree", path+"PKUTree_gKK-3000_R-1500_"+year+".root");

            if variable not in [str(i) for i in df_QCD.GetColumnNames()]:
                Intime_plot_var = "Var1"
                try:
                    print "\n"
                    print eval(variable)
                    if MODE in ["MCvsDATA","COMP"]:
                        df_data= df_data.Define(   "Var1",eval(variable))
                    df_QCD     = df_QCD.Define(    "Var1",eval(variable))
                    df_WJets   = df_WJets.Define(  "Var1",eval(variable))
                    df_TTbar   = df_TTbar.Define(  "Var1",eval(variable))
                    df_STop    = df_STop.Define(   "Var1",eval(variable))
                    df_Rest    = df_Rest.Define(   "Var1",eval(variable))
                    df_Signal1 = df_Signal1.Define("Var1",eval(variable))
                    df_Signal2 = df_Signal2.Define("Var1",eval(variable))
                except (NameError,SyntaxError):
                    print "\n"
                    print variable
                    if MODE in ["MCvsDATA","COMP"]:
                        df_data    = df_data.Define(   "Var1","return (%s);"%(variable))
                    df_QCD     = df_QCD.Define(    "Var1","return (%s);"%(variable))
                    df_WJets   = df_WJets.Define(  "Var1","return (%s);"%(variable))
                    df_TTbar   = df_TTbar.Define(  "Var1","return (%s);"%(variable))
                    df_STop    = df_STop.Define(   "Var1","return (%s);"%(variable))
                    df_Rest    = df_Rest.Define(   "Var1","return (%s);"%(variable))
                    df_Signal1 = df_Signal1.Define("Var1","return (%s);"%(variable))
                    df_Signal2 = df_Signal2.Define("Var1","return (%s);"%(variable))
            else:
                Intime_plot_var = variable

            for Added_var in self.Intime_Cut_Variable:
                if Added_var in cut:
                    if Added_var not in [str(i) for i in df_QCD.GetColumnNames()]:
                        if MODE in ["MCvsDATA","COMP"]:
                            df_data= df_data.Define(Added_var , eval(Added_var))
                        df_QCD     = df_QCD.Define(Added_var  , eval(Added_var))
                        df_WJets   = df_WJets.Define(Added_var, eval(Added_var))
                        df_TTbar   = df_TTbar.Define(Added_var, eval(Added_var))
                        df_STop    = df_STop.Define(Added_var , eval(Added_var))
                        df_Rest    = df_Rest.Define(Added_var , eval(Added_var))
                        df_Signal1 = df_Signal1.Define(Added_var, eval(Added_var))
                        df_Signal2 = df_Signal2.Define(Added_var, eval(Added_var))


        if MODE in ["DECO"]:
            for c in ["R4q","R3q","R2q","W","Rlqq","Rlq","g",  "Rest", "Neut", "qlea"]:
                exec('h_Signal1_%s = TH1D("h_Signal1_%s","h_Signal1_%s"'%(c,c,c)+'";%s;%s"'%(xtitle,ytitle)+',nbin,min,max); h_Signal1_%s.Sumw2();'%(c));

        if MODE in ["COMP"]:
            print " ---> define histograms for SR, CR?";
            h_DMR    = TH1D("h_DMR"    ,"h_DMR"  +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_DMR.Sumw2();
            h_QCD    = TH1D("h_QCD"    ,"h_QCD"  +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_QCD.Sumw2();
            h_QCD_S  = TH1D("h_QCD_S"  ,"h_QCD_S"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_QCD_S.Sumw2();
            h_TTbar  = TH1D("h_TTbar"  ,"h_TTbar"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_TTbar.Sumw2();
            h_STop   = TH1D("h_STop"   ,"h_STop" +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_STop.Sumw2();
            h_Rest   = TH1D("h_Rest"   ,"h_Rest" +";%s;%s"%(xtitle,ytitle),nbin,min,max); h_Rest.Sumw2();
            h_TotalMC= TH1D("h_TotalMC","h_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_TotalMC.Sumw2();
            #c_data   = h_data.Clone("c_data");
            #c_Signal1= h_Signal1.Clone("c_Signal1");  
            #c_Signal2= h_Signal2.Clone("c_Signal2");  #c_Signal3=h_Signal3.Clone("c_Signal3");  c_Signal4=h_Signal4.Clone("c_Signal4");
            #c_QCD    = h_QCD.Clone("c_QCD");  c_WJets=h_WJets.Clone("c_WJets");  c_TTbar=h_TTbar.Clone("c_TTbar");  c_STop=h_STop.Clone("c_STop");  c_Rest=h_Rest.Clone("c_Rest"); c_TotalMC=h_TotalMC.Clone("c_TotalMC");

        hstack_TotalMC= THStack("hstack_TotalMC","hstack_TotalMC"+";%s;%s"%(xtitle,ytitle));                


        #=================== SET WEIGHTS, SCALE TREES, DEFINE TOTAL AND STACK  =================================================
        weight="weight";
        if SFs==1: weight="weight_center";

        if MODE=="DECO":self.Signal_Scale1=1;

        if Intime:
            Bin = ("Plot",";%s;%s"%(xtitle,ytitle),nbin,min,max)
            if MODE in ["MCvsDATA","COMP"]:
                h_data = df_data.Filter(  cut).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_data" );
            h_QCD     = df_QCD.Filter(    cut).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_QCD"  );
            h_WJets   = df_WJets.Filter(  cut).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_WJets");
            h_TTbar   = df_TTbar.Filter(  cut).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_TTbar");
            h_STop    = df_STop.Filter(   cut).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_STop" );
            h_Rest    = df_Rest.Filter(   cut).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Rest" );
            h_Signal1 = df_Signal1.Filter(cut).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1"); h_Signal1.Scale(self.Signal_Scale1);
            h_Signal2 = df_Signal2.Filter(cut).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal2"); h_Signal2.Scale(self.Signal_Scale2);
            h_QCD_S   = df_QCD.Filter(   cut1).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_QCD_S"  );
        else:
            for region in ["h_"]:
                for sample in ["data","Signal1","Signal2","QCD","WJets","TTbar","STop","Rest"]:
                    if MODE in ["MCvsDATA","COMP"]:
                        if sample in ["data"]                         : eval("t_"+sample).Draw("(%s) >> %s%s"%(variable,region,sample)   ,                                           cut  );   # No weights on data
                    eval(                                                    "t_Signal1").Draw("(%s) >> %s%s"%(variable,region,"Signal1"), "(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut) );   # Extra scaling on signal
                    eval(                                                    "t_Signal2").Draw("(%s) >> %s%s"%(variable,region,"Signal2"), "(%s*%s)*(%s)"%(weight,self.Signal_Scale2,cut) );   # Extra scaling on signal
                    if sample in ["QCD","WJets","TTbar","STop","Rest"]: eval("t_"+sample).Draw("(%s) >> %s%s"%(variable,region,sample)   , "(%s*%s)*(%s)"%(weight,1                 ,cut) );   # MC only weights on MC BKG
                if MODE=="COMP":
                    eval("t_QCD").Draw("(%s) >> %s%s"%(variable,region,"QCD_S"), "(%s*%s)*(%s)"%(weight,1  ,cut1) );  print "  got SR QCD histo"; 

        if MODE in ["DECO"]:  #  comp = ["R4q","R3q","R2q","W","Rlqq","Rlq","Rest"];
            self.Signal_Scale1=1;
            R4q  = "( R4q_a   == 1)";
            R3q  = "( R3q_a   == 1)";
            R2q  = "( R2q_a   == 1)";
            W    = "( w_a     == 1)";
            Rlqq = "( Rlqq_a  == 1)";
            Rlq  = "( Rlq_a   == 1)";
            g    = "( gKK_g_a == 1)";
            Neut = "( R3q_taudecay_a == 1 || Rlqq_a == 1 || Rlq_a == 1 )";
            Rest = "(!"+R4q+" && !"+R3q+" && !"+R2q+" && !"+W+" && !"+Rlqq+" && !"+Rlq+" && !"+g+" )";
            qlea = "( !"+Rest+" && !"+Neut+" &&  ("+R3q+" || "+R2q+" || "+W+" ))";

            #Wa_MO_AK8 == 1 :  Max score is the max mass jet
            #Wa_MO_AK8 == 2 :  Max score is the mid mass jet
            #Wa_MO_AK8 == 3 :  Max score is the min mass jet
            if Intime:
                h_Signal1_R4q  = df_Signal1.Filter( cut+"&&" +R4q ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_R4q" ); h_Signal1_R4q.Scale(self.Signal_Scale1);
                h_Signal1_R3q  = df_Signal1.Filter( cut+"&&" +R3q ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_R3q" ); h_Signal1_R3q.Scale(self.Signal_Scale1);
                h_Signal1_R2q  = df_Signal1.Filter( cut+"&&" +R2q ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_R2q" ); h_Signal1_R2q.Scale(self.Signal_Scale1);
                h_Signal1_W    = df_Signal1.Filter( cut+"&&"   +W ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_W"   ); h_Signal1_W.Scale(self.Signal_Scale1);
                h_Signal1_Rlqq = df_Signal1.Filter( cut+"&&"+Rlqq ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_Rlqq"); h_Signal1_Rlqq.Scale(self.Signal_Scale1);
                h_Signal1_Rlq  = df_Signal1.Filter( cut+"&&" +Rlq ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_Rlq" ); h_Signal1_Rlq.Scale(self.Signal_Scale1);
                h_Signal1_g    = df_Signal1.Filter( cut+"&&"   +g ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_g"   ); h_Signal1_g.Scale(self.Signal_Scale1);
                h_Signal1_Rest = df_Signal1.Filter( cut+"&&"+Rest ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_Rest"); h_Signal1_Rest.Scale(self.Signal_Scale1);
                h_Signal1_Neut = df_Signal1.Filter( cut+"&&"+Neut ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_Neut"); h_Signal1_Neut.Scale(self.Signal_Scale1);
                h_Signal1_qlea = df_Signal1.Filter( cut+"&&"+qlea ).Histo1D(Bin, Intime_plot_var, weight).GetValue().Clone("h_Signal1_qlea"); h_Signal1_qlea.Scale(self.Signal_Scale1);
            else:
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_R4q") ,"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+R4q ) );
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_R3q") ,"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+R3q ) );
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_R2q") ,"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+R2q ) );
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_W")   ,"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+W   ) );
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_Rlqq"),"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+Rlqq) );
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_Rlq" ),"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+Rlq ) );
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_g")   ,"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+g   ) );
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_Rest"),"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+Rest) );
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_Neut"),"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+Neut) );
                eval("t_Signal1").Draw("(%s) >> %s"%(variable,"h_Signal1_qlea"),"(%s*%s)*(%s)"%(weight,self.Signal_Scale1,cut+"&&"+qlea) );

            h_TotalS1=TH1D("h_TotalS1","h_TotalS1"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_TotalS1.Sumw2();
            for h in [           h_Signal1_R4q, h_Signal1_R3q, h_Signal1_R2q, h_Signal1_W, h_Signal1_Rlqq, h_Signal1_Rlq, h_Signal1_g,   h_Signal1_Rest                                 ]: h_TotalS1.Add(h);
            for h in [h_TotalS1, h_Signal1_R4q, h_Signal1_R3q, h_Signal1_R2q, h_Signal1_W, h_Signal1_Rlqq, h_Signal1_Rlq, h_Signal1_g,   h_Signal1_Rest,  h_Signal1_Neut, h_Signal1_qlea]: h=UnderOverFlow1D(h);

        hstack_TotalMC.Add(h_Rest);hstack_TotalMC.Add(h_STop);hstack_TotalMC.Add(h_TTbar);hstack_TotalMC.Add(h_WJets);hstack_TotalMC.Add(h_QCD);
        h_TotalMC = TH1D("h_TotalMC","h_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max); h_TotalMC.Sumw2(); h_TotalMC.Add(h_QCD); h_TotalMC.Add(h_WJets); h_TotalMC.Add(h_TTbar); h_TotalMC.Add(h_STop); h_TotalMC.Add(h_Rest);
        h_data.SetBinErrorOption(TH1D.kPoisson);

        for h in [h_data,h_Signal1,h_Signal2,h_QCD,h_WJets,h_TTbar,h_STop,h_Rest,h_TotalMC ]:h=UnderOverFlow1D(h);

        # NORMALIZE SIGNAL TO BKG 
        NORM_s1=h_TotalMC.Integral()/(h_Signal1.Integral()+0.000001);
        if REGION[:3] in ["PS3","CR4","SR4"]                        : NORM_s1=0;
        elif                     NORM_s1>=300000: NORM_s1=100000;
        elif  300000>NORM_s1 and NORM_s1>=100000: NORM_s1= 30000;
        elif  100000>NORM_s1 and NORM_s1>= 30000: NORM_s1= 10000;
        elif   30000>NORM_s1 and NORM_s1>= 10000: NORM_s1=  3000;
        elif   10000>NORM_s1 and NORM_s1>=  3000: NORM_s1=  1000;
        elif    3000>NORM_s1 and NORM_s1>=  1000: NORM_s1=   300;
        elif    1000>NORM_s1 and NORM_s1>=   300: NORM_s1=   100;
        elif     300>NORM_s1 and NORM_s1>=   100: NORM_s1=    30;
        elif     100>NORM_s1 and NORM_s1>=    30: NORM_s1=    10;
        elif      30>NORM_s1 and NORM_s1>=    10: NORM_s1=     3;
        else                                    : NORM_s1=     1;
        NORM_s2=h_TotalMC.Integral()/(h_Signal2.Integral()+0.000001);
        if REGION[:3] in ["PS2","CR1","CR2","CR3","SR1","SR2","SR3"]: NORM_s2=0;
        if                       NORM_s2>=300000: NORM_s2=     0; 
        elif                     NORM_s2>=300000: NORM_s2=100000;
        elif  300000>NORM_s2 and NORM_s2>=100000: NORM_s2= 30000;
        elif  100000>NORM_s2 and NORM_s2>= 30000: NORM_s2= 10000;
        elif   30000>NORM_s2 and NORM_s2>= 10000: NORM_s2=  3000;
        elif   10000>NORM_s2 and NORM_s2>=  3000: NORM_s2=  1000;
        elif    3000>NORM_s2 and NORM_s2>=  1000: NORM_s2=   300;
        elif    1000>NORM_s2 and NORM_s2>=   300: NORM_s2=   100;
        elif     300>NORM_s2 and NORM_s2>=   100: NORM_s2=    30;
        elif     100>NORM_s2 and NORM_s2>=    30: NORM_s2=    10;
        elif      30>NORM_s2 and NORM_s2>=    10: NORM_s2=     3;
        else                                    : NORM_s2=     1;
        if MODE=="MC":
            h_Signal1.Scale(NORM_s1);
            h_Signal2.Scale(NORM_s2);


        #     Define Significance
        #-------- build denominator of significance ------------
        h_SqrtTotalMC=h_TotalMC.Clone("h_SqrtTotalMC"); h_SqrtTotalMC.Sumw2();

        for i in range(1,h_SqrtTotalMC.GetNbinsX()+1,1): 
            h_SqrtTotalMC.SetBinContent(i, TMath.Sqrt( h_SqrtTotalMC.GetBinContent(i)+1 ) );
            h_SqrtTotalMC.SetBinError( i, h_TotalMC.GetBinError(i)/(2*TMath.Sqrt(h_SqrtTotalMC.GetBinContent(i)+1) ) );
        h_Signif1=h_Signal1.Clone("h_Signif1"); h_Signif1.Divide( h_SqrtTotalMC );
        h_SoB1   =h_Signal1.Clone("h_SoB1");    h_SoB1.Divide( h_TotalMC );
        h_Signif2=h_Signal2.Clone("h_Signif2"); h_Signif2.Divide( h_SqrtTotalMC );

        # Significance Ordered Transformation
        #for h in [h_data,h_Signal1,h_Signal2,h_QCD,h_WJets,h_TTbar,h_STop,h_Rest,h_TotalMC ] : h=SOT_1D(h,h_Signif1);


        #if MODE in ["COMP"]:
        #    c_TotalMC = TH1D("c_TotalMC","c_TotalMC"+";%s;%s"%(xtitle,ytitle),nbin,min,max);c_TotalMC.Sumw2();
        #    c_TotalMC.Add(c_QCD);c_TotalMC.Add(c_WJets);c_TotalMC.Add(c_TTbar);c_TotalMC.Add(c_STop);c_TotalMC.Add(c_Rest);
        #    c_REST    = TH1D("c_REST","c_REST"+";%s;%s"%(xtitle,ytitle),nbin,min,max);c_REST.Sumw2(); c_REST.Add(c_WJets);c_REST.Add(c_TTbar);c_REST.Add(c_STop);c_REST.Add(c_Rest); ## that is for non QCD in CR
        #    c_data.SetBinErrorOption(TH1D.kPoisson);
        #    for h in [ c_data,c_Signal1,c_Signal2,c_QCD,c_WJets,c_TTbar,c_STop,c_TotalMC, c_Rest ]:h=UnderOverFlow1D(h);   # UNDEROVERFLOWS for All histos used after reweighting

            #h_TotalMC.Scale(norm); h_QCD.Scale(norm);  h_WJets.Scale(norm);  h_TTbar.Scale(norm);  h_STop.Scale(norm);  h_Rest.Scale(norm);

        #if options.FBT:  # Transformes the spectra in such a way where BKG is flat!
        #    for h in [h_data,h_Signal1,h_Signal2,h_QCD,h_WJets,h_TTbar,h_STop,h_Rest,h_TotalMC ]: h=FBT_(h_TotalMC, h);

        norm=h_data.Integral()/(h_TotalMC.Integral()+0.00001); 
        #print "  norm=",norm;

        if MODE in ["MC","MCvsDATA"]: #---------- histogram cosmetics ---------------------
            h_data.SetLineColor(self.color_palet["data"]); h_data.SetFillColor(self.color_palet["data"]);
            h_Signal1.SetLineColor(4);       h_Signal1.SetFillStyle(0); h_Signal1.SetLineWidth(3); h_Signal1.SetLineStyle(7);
            h_Signif1.SetLineColor(4);       h_Signif1.SetFillStyle(0); h_Signif1.SetLineWidth(3); h_Signif1.SetLineStyle(7);
            h_SoB1.SetLineColor(7);          h_SoB1.SetFillStyle(0);    h_SoB1.SetLineWidth(3);    h_SoB1.SetLineStyle(3);
            h_Signal2.SetLineColor(kGreen+2);h_Signal2.SetFillStyle(0); h_Signal2.SetLineWidth(3); h_Signal2.SetLineStyle(7);
            h_Signif2.SetLineColor(kGreen+2);h_Signif2.SetFillStyle(0); h_Signif2.SetLineWidth(3); h_Signif2.SetLineStyle(7);
            #h_Signal3.SetLineColor(1);       h_Signal3.SetFillStyle(0); h_Signal3.SetLineWidth(3); h_Signal3.SetLineStyle(7);
            #h_Signal4.SetLineColor(93);      h_Signal4.SetFillStyle(0); h_Signal4.SetLineWidth(3); h_Signal4.SetLineStyle(7);
            h_QCD.SetLineColor(0);     h_QCD.SetFillColor(  self.color_palet["QCD"]  ); h_QCD.SetLineWidth(0);
            h_WJets.SetLineColor(0);   h_WJets.SetFillColor(self.color_palet["WJets"]); h_WJets.SetLineWidth(0);
            h_TTbar.SetLineColor(0);   h_TTbar.SetFillColor(self.color_palet["TTbar"]); h_TTbar.SetLineWidth(0);
            h_STop.SetLineColor(0);    h_STop.SetFillColor( self.color_palet["STop"] ); h_STop.SetLineWidth(0);
            h_Rest.SetLineColor(0);    h_Rest.SetFillColor( self.color_palet["Rest"]);  h_Rest.SetLineWidth(0);
            h_TotalMC.SetLineStyle(3); h_TotalMC.SetMarkerStyle(0); h_TotalMC.SetLineWidth(5); h_TotalMC.SetLineColor(15);
        if MODE in ["DECO"]:
            h_TotalMC.SetLineStyle(1);      h_TotalMC.SetMarkerStyle(0);   h_TotalMC.SetLineWidth(3);     h_TotalMC.SetLineColor(1);
            h_Signal1.SetLineColor(4);      h_Signal1.SetFillStyle(0);     h_Signal1.SetLineWidth(3);     h_Signal1.SetLineStyle(7);
            h_TotalS1.SetLineColor(2);      h_TotalS1.SetFillStyle(0);     h_TotalS1.SetLineWidth(3);     h_TotalS1.SetLineStyle(3);
            h_Signal1_R4q.SetLineColor(8);  h_Signal1_R4q.SetFillStyle(0); h_Signal1_R4q.SetLineWidth(3); h_Signal1_R4q.SetLineStyle(7);
            h_Signal1_R3q.SetLineColor(92); h_Signal1_R3q.SetFillStyle(0); h_Signal1_R3q.SetLineWidth(3); h_Signal1_R3q.SetLineStyle(7);
            h_Signal1_R2q.SetLineColor(46); h_Signal1_R2q.SetFillStyle(0); h_Signal1_R2q.SetLineWidth(3); h_Signal1_R2q.SetLineStyle(7);
            h_Signal1_W.SetLineColor(4);    h_Signal1_W.SetFillStyle(0);   h_Signal1_W.SetLineWidth(3);   h_Signal1_W.SetLineStyle(7);
            h_Signal1_Rlqq.SetLineColor(65);h_Signal1_Rlqq.SetFillStyle(0);h_Signal1_Rlqq.SetLineWidth(3);h_Signal1_Rlqq.SetLineStyle(7);
            h_Signal1_Rlq.SetLineColor(15); h_Signal1_Rlq.SetFillStyle(0); h_Signal1_Rlq.SetLineWidth(3); h_Signal1_Rlq.SetLineStyle(7);
            h_Signal1_g.SetLineColor(96);   h_Signal1_g.SetFillStyle(0);   h_Signal1_g.SetLineWidth(3);   h_Signal1_g.SetLineStyle(7);
            h_Signal1_Rest.SetLineColor(18);h_Signal1_Rest.SetFillStyle(0);h_Signal1_Rest.SetLineWidth(3);h_Signal1_Rest.SetLineStyle(7);
            h_Signal1_Neut.SetLineColor(64);h_Signal1_Neut.SetFillStyle(0);h_Signal1_Neut.SetLineWidth(3);h_Signal1_Neut.SetLineStyle(7);
            h_Signal1_qlea.SetLineColor(95);h_Signal1_qlea.SetFillStyle(0);h_Signal1_qlea.SetLineWidth(3);h_Signal1_qlea.SetLineStyle(3);
        if MODE=="COMP":
            print " ---> Forming DMR, scaling histos & cosmetics";
            h_DMR=h_data.Clone("h_data");  h_DMR.Add(h_WJets,-1);h_DMR.Add(h_TTbar,-1);h_DMR.Add(h_Rest,-1);h_DMR.Add(h_STop,-1);    h_QCD_S=UnderOverFlow1D(h_QCD_S); 
            #for h in [h_DMR,h_QCD_S,h_data,h_Signal1,h_Signal2,h_QCD,h_WJets,h_TTbar,h_STop,h_Rest,h_TotalMC]:h=UnderOverFlow1D(h);   # UNDEROVERFLOWS for All histos used after reweighting

            Norm =h_QCD_S.Integral()/h_QCD.Integral(); print "Norm = h_QCD_S.Integral()/h_QCD.Integral()=",h_QCD_S.Integral(),"/",h_QCD.Integral();
            Norm2=h_QCD_S.Integral()/h_DMR.Integral(); print "Norm2= h_QCD_S.Integral()/h_DMR.Integral()=",h_QCD_S.Integral(),"/",h_DMR.Integral();
            NORM =h_DMR.Integral()/h_QCD.Integral();   print "NORM = h_DMR.Integral()/h_QCD.Integral();=" ,h_DMR.Integral()  ,"/",h_QCD.Integral();
            h_QCD_S.Scale(1*NORM);
            h_QCD.Scale(Norm*NORM);
            h_DMR.Scale(Norm2*NORM);
            h_DMR.SetLineColor(1);   h_DMR.SetFillStyle(0); h_DMR.SetLineWidth(3); h_DMR.SetLineStyle(0);             h_DMR.SetMarkerStyle(1); h_DMR.SetMarkerSize(3);
            h_QCD_S.SetLineColor(2);   h_QCD_S.SetFillStyle(0); h_QCD_S.SetLineWidth(3); h_QCD_S.SetLineStyle(0);     h_QCD_S.SetMarkerStyle(1); h_QCD_S.SetMarkerSize(7);
            h_QCD.SetLineColor(kBlue);   h_QCD.SetFillStyle(0); h_QCD.SetLineWidth(3); h_QCD.SetLineStyle(0);         h_QCD.SetMarkerStyle(0); h_QCD.SetMarkerSize(0);



        #============ DRAW TOP PAD =====================
        if MODE in ["MC"]:
            h_TotalMC.Draw("e"); h_TotalMC.GetXaxis().SetNdivisions(509);
            hstack_TotalMC.Draw("same HIST"); # For unc-bars
            h_TotalMC.Draw("same e"   ); #print "last bin",h_TotalMC.GetBinContent(27);
            if REGION[:3] in ["PS2","CR1","CR2","CR3","SR1","SR2","SR3"]:h_Signal1.Draw("same HIST");
            if REGION[:3] in ["PS3","CR4","SR4"]                        :h_Signal2.Draw("same HIST");
            #h_Signal3.Draw("same HIST");
            #h_Signal4.Draw("same HIST");
            #canvas_controlplot.Update(); 

        if MODE in ["MCvsDATA"]:
            h_TotalMC.Draw("e"); h_TotalMC.GetXaxis().SetNdivisions(509);   
            h_data.Draw("e same"); h_data.GetXaxis().SetNdivisions(509);
            hstack_TotalMC.Draw("same HIST"); # For unc-bars
            h_data.Draw("same e");  #!needed 2nd time to draw data
            h_TotalMC.Draw("same e"   );
            if REGION[:3] in ["PS2","CR1","CR2","CR3","SR1","SR2","SR3"]:h_Signal1.Draw("same HIST");
            if REGION[:3] in ["PS3","CR4","SR4"]                        :h_Signal2.Draw("same HIST");
            #h_Signal3.Draw("same HIST");
            #h_Signal4.Draw("same HIST");
            canvas_controlplot.Update(); 

        if MODE in ["DECO"]:
            h_Signal1.Scale( 1/self.Signal_Scale1);
            h_TotalMC.Scale( h_Signal1.Integral()/(h_TotalMC.Integral()+0.000001) );
            if logy==0:h_Signal1.GetYaxis().SetRangeUser(  0,1.3*TMath.Max(h_Signal1.GetMaximum(),h_TotalMC.GetMaximum())  );
            if logy==1:h_Signal1.GetYaxis().SetRangeUser(0.1, 30*h_Signal1.GetMaximum());
            h_Signal1.Draw("HIST");
            h_TotalMC.Draw("same HIST");
            ###h_TotalS1.Draw("same HIST");
            #h_Signal1_R4q.Draw("same HIST");
            #h_Signal1_R3q.Draw("same HIST");
            #h_Signal1_R2q.Draw("same HIST");
            #h_Signal1_W.Draw("same HIST");
            #h_Signal1_Rlqq.Draw("same HIST");
            #h_Signal1_Rlq.Draw("same HIST");
            #h_Signal1_Rest.Draw("same HIST");
            #h_Signal1_g.Draw("same HIST");
            h_Signal1_Neut.Draw("same HIST");
            #h_Signal1_qlea.Draw("same HIST");
            if logy==1:canvas_controlplot.SetLogy();
            canvas_controlplot.Update(); 

        if MODE in ["COMP"]:
            print " ---> Drawing Histos";
            h_QCD_S.Draw("HIST,e");     print "h_QCD_S",h_QCD_S.Integral();
            h_QCD.Draw("same HIST,e");  print "h_QCD"  ,h_QCD.Integral();
            h_DMR.Draw("same HIST,e");  print "h_DMR"  ,h_DMR.Integral();

        #if logy and MODE!="DECO":
        #    fPads1.SetLogy();
        #    fPads1.Update();
        #    canvas_controlplot.Update();
            #if MODE is "MC":
            #    fPads2.SetLogy();
            #    fPads2.Update();
            #    canvas_controlplot.Update(); 


        #---------------- Add text in top pad -----------------------
        banner          = TLatex(0.96,0.96,str(lumi)+" fb^{-1} (13 TeV)");   banner.SetNDC();   banner.SetTextSize(0.034);     banner.SetTextFont(42);    banner.SetTextAlign(31);    banner.SetLineWidth(2);    banner.Draw();
        CMS             = TLatex(0.22,0.96,"CMS"                        );      CMS.SetNDC();      CMS.SetTextSize(0.042);        CMS.SetTextFont(42);       CMS.SetTextAlign(31);       CMS.SetLineWidth(2);       CMS.Draw();
        if MODE=="MCvsData":
            Extratext   = TLatex(0.24,0.96,"Preliminary"                );Extratext.SetNDC();Extratext.SetTextSize(0.034);  Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        if MODE=="DECO" or MODE=="MC":
            Extratext   = TLatex(0.24,0.96,"Simulation"                 );Extratext.SetNDC();Extratext.SetTextSize(0.034);  Extratext.SetTextFont(52); Extratext.SetTextAlign(11); Extratext.SetLineWidth(2); Extratext.Draw();
        RegionTxt       = TLatex(0.15,0.88,"%s"%(REGION)                );RegionTxt.SetNDC();RegionTxt.SetTextSize(0.042);  RegionTxt.SetTextFont(42);    RegionTxt.SetLineWidth(2); RegionTxt.Draw();
        if MODE in ["MCvsDATA","COMP"]:
            D_o_MC_txt  = TLatex(0.15,0.83,"Data / MC = %.2f"%(norm)    );D_o_MC_txt.SetNDC();D_o_MC_txt.SetTextSize(0.042);D_o_MC_txt.SetTextFont(42);   D_o_MC_txt.SetLineWidth(2);D_o_MC_txt.Draw();
        if SFs :
            SFsCorr     = TLatex(0.55, 0.96, "SFs corrected"            );   SFsCorr.SetNDC();   SFsCorr.SetTextSize(0.034);   SFsCorr.SetTextFont(52);   SFsCorr.SetTextAlign(11);   SFsCorr.SetLineWidth(2);   SFsCorr.Draw();
        #canvas_controlplot.Update(); 


        #========== DRAW BOTTOM PAD ============================================
        if MODE in ["MCvsDATA"] : #--------- Data / MC on 2nd pad ---------------------
            fPads2.cd(); 
            h_Ratio = h_data.Clone("h_Ratio"); h_Ratio.Divide( h_TotalMC ); MaxY=2; #TMath.Max( 2,  TMath.Min(3,h_Ratio.GetMaximum()*1.1) );
            h_Ratio.SetLineColor(1); h_Ratio.SetLineWidth(2); h_Ratio.SetMarkerStyle(8); h_Ratio.SetMarkerSize(0.7); h_Ratio.GetYaxis().SetRangeUser( 0 , MaxY );  h_Ratio.GetYaxis().SetNdivisions(504,0);
            h_Ratio.GetYaxis().SetTitle("Data / MC  ");  h_Ratio.GetYaxis().SetTitleOffset(0.35);  h_Ratio.GetYaxis().SetTitleSize(0.13);  h_Ratio.GetYaxis().SetTitleSize(0.13);  h_Ratio.GetYaxis().SetLabelSize(0.11); h_Ratio.GetXaxis().SetLabelSize(0.1); h_Ratio.GetXaxis().SetTitleOffset(0.7); h_Ratio.GetXaxis().SetTitleSize(0.14); 
            axis1=TGaxis( min,1,max,1, 0,0,0, "L"); axis1.SetLineColor(1); axis1.SetLineWidth(1);  #axis1->SetLabelColor(16); #fPads2.SetGridx(); #fPads2.SetGridy();
            h_Ratio=RationUnc(h_data,h_TotalMC,h_Ratio,MaxY);
            h_Ratio.Draw("e0"); axis1.Draw();
            fPads2.RedrawAxis(); fPads2.Update();
            fPads1.RedrawAxis(); fPads1.Update();

        #--------- Significances on 2nd pad ---------------------
        if MODE in ["MC"]:  
            fPads2.cd(); MaxY=7; MinY=0; 
            #fPads2.SetLogy();MinY=0.01;
            gStyle.SetPadTickY(1);
            #axis2=TGaxis( min,2,max,2, 0,0,0, "L"); axis2.SetLineColor(2); axis2.SetLineWidth(1);axis2.Draw();
            #axis3=TGaxis( min,5,max,5, 0,0,0, "L"); axis3.SetLineColor(3); axis3.SetLineWidth(1);axis3.Draw();

            h_SqrtTotalMC=h_TotalMC.Clone("h_SqrtTotalMC"); h_SqrtTotalMC.Sumw2();
            for i in range(1,h_SqrtTotalMC.GetNbinsX()+1,1): h_SqrtTotalMC.SetBinContent(i, TMath.Sqrt( 1.0*h_SqrtTotalMC.GetBinContent(i)+1 ) ); h_SqrtTotalMC.SetBinError( i, TMath.Sqrt( 1.0*h_SqrtTotalMC.GetBinError(i) ) );
            h_Signif1=h_Signal1.Clone("h_Signif1"); h_Signif1.Divide( h_SqrtTotalMC );  h_SoB1=h_Signal1.Clone("h_SoB1"); h_SoB1.Divide( h_TotalMC ); #h_SqrtTotalMC
            h_Signif2=h_Signal2.Clone("h_Signif2"); h_Signif2.Divide( h_SqrtTotalMC );  h_SoB2=h_Signal2.Clone("h_SoB2"); h_SoB2.Divide( h_TotalMC ); #h_SqrtTotalMC

            h_Signif1.GetYaxis().SetNdivisions(504,1);
            h_Signif2.GetYaxis().SetNdivisions(504,1);

            h_Signif1.GetYaxis().SetTitle(" S / #sqrt{B+1} "); h_Signif1.GetYaxis().SetTitleOffset(0.38); h_Signif1.GetYaxis().SetTitleSize(0.13); h_Signif1.GetYaxis().SetLabelSize(0.11); h_Signif1.GetXaxis().SetLabelSize(0.1);h_Signif1.GetXaxis().SetTitleOffset(0.7);h_Signif1.GetXaxis().SetTitleSize(0.14);

            MaxY = TMath.Max(MaxY, h_Signif1.GetMaximum()); MaxY = TMath.Max(MaxY, h_Signif2.GetMaximum()); #MaxY=TMath.Max(MaxY,h_Signif3.GetMaximum());  MaxY=TMath.Max(MaxY,h_Signif4.GetMaximum());  
            MaxY = TMath.Max(MaxY, h_SoB1.GetMaximum());
            MaxY=MaxY*1.3;
            h_Signif1.GetYaxis().SetRangeUser(MinY,MaxY); h_SoB1.GetYaxis().SetRangeUser(MinY,MaxY);
            h_Signif1.Draw("hist");  #h_SoB1.Scale( MaxY ); h_SoB1.SetLineColor(7);h_SoB1.SetLineStyle(3);  #h_SoB1.Draw("same hist");
            h_Signif2.Draw("hist,same"); #h_Signif3.Draw("hist,same");h_Signif4.Draw("hist,same"); 
            #axis2.Draw();axis3.Draw();
            #ResultList = [ LeftCutBin , RightCutBin , S.GetBinLowEdge(LeftCutBin) , S.GetBinLowEdge(RightCutBin+1), str(BKGrjc), str(Sig_Eff), str(SigMax_Print) ];  
            if self.Optimal:
                length=0.06*(h_TotalMC.GetBinLowEdge(h_TotalMC.GetNbinsX()+1)-h_TotalMC.GetBinLowEdge(1)); 
                if REGION[:3] in ["PS2","CR1","CR2","CR3","SR1","SR2","SR3"]:
                    V1 = OptimalCut(h_TotalMC, h_Signal1, NORM_s1); 
                    sig1_Leftcut =TGaxis( V1[2], MinY, V1[2],MaxY, 0,0,0, "L"); sig1_Leftcut.SetLineColor(4);   sig1_Leftcut.SetLineWidth(3);  sig1_Leftcut.Draw();
                    #sig1_Midcut  =TGaxis( V1[7], MinY, V1[7],MaxY, 0,0,0, "L"); sig1_Midcut.SetLineColor(4);    sig1_Midcut.SetLineWidth(3);   sig1_Midcut.Draw();
                    sig1_Rightcut=TGaxis( V1[3], MinY, V1[3],MaxY, 0,0,0, "L"); sig1_Rightcut.SetLineColor(4);  sig1_Rightcut.SetLineWidth(3); sig1_Rightcut.Draw();
                    PrintValus1 = TLatex( h_TotalMC.GetBinLowEdge(2), MaxY*0.75,"    r "+V1[4]+"%   e "+V1[5]+"%   si "+V1[6]+"%   S/srB "+V1[7]); PrintValus1.SetTextColor(4); PrintValus1.SetTextSize(0.13); PrintValus1.Draw();
                    ar_L = TArrow( V1[2] ,MaxY/3, V1[2]+length, MaxY/3, 0.03 , "|>"); ar_L.SetLineWidth(3); ar_L.SetFillColor(4); ar_L.SetLineColor(4); ar_L.SetAngle(30); ar_L.Draw();
                    ar_R = TArrow( V1[3] ,MaxY/3, V1[3]-length, MaxY/3, 0.03 , "|>"); ar_R.SetLineWidth(3); ar_R.SetFillColor(4); ar_R.SetLineColor(4); ar_R.SetAngle(30); ar_R.Draw();
                    if abs(V1[2]-0.5)<abs(V1[3]-0.5): self.exponent = V1[2]; 
                    if abs(V1[2]-0.5)>abs(V1[3]-0.5): self.exponent = V1[3];
                if REGION[:3] in ["PS3","CR4","SR4"]:
                    V2 = OptimalCut(h_TotalMC, h_Signal2, NORM_s2);
                    sig2_Leftcut =TGaxis( V2[2], MinY, V2[2],MaxY, 0,0,0, "L"); sig2_Leftcut.SetLineColor(kGreen+2);   sig2_Leftcut.SetLineWidth(3);  sig2_Leftcut.Draw();
                    #sig2_Midcut  =TGaxis( V2[7], MinY, V2[7],MaxY, 0,0,0, "L"); sig2_Midcut.SetLineColor(kGreen+2);    sig2_Midcut.SetLineWidth(3);   sig2_Midcut.Draw();
                    sig2_Rightcut=TGaxis( V2[3], MinY, V2[3],MaxY, 0,0,0, "L"); sig2_Rightcut.SetLineColor(kGreen+2);  sig2_Rightcut.SetLineWidth(3); sig2_Rightcut.Draw();
                    PrintValus2 = TLatex( h_TotalMC.GetBinLowEdge(2), MaxY*0.55,"    r "+V2[4]+"%   e "+V2[5]+"%   si "+V2[6]+"%   S/srB "+V2[7]); PrintValus2.SetTextColor(kGreen+2); PrintValus2.SetTextSize(0.13); PrintValus2.Draw();
                    ar_l = TArrow( V2[2] ,MaxY/4, V2[2]+length, MaxY/4, 0.03 , "|>"); ar_l.SetLineWidth(3); ar_l.SetFillColor(kGreen+2); ar_l.SetLineColor(kGreen+2); ar_l.SetAngle(30); ar_l.Draw();
                    ar_r = TArrow( V2[3] ,MaxY/4, V2[3]-length, MaxY/4, 0.03 , "|>"); ar_r.SetLineWidth(3); ar_r.SetFillColor(kGreen+2); ar_r.SetLineColor(kGreen+2); ar_r.SetAngle(30); ar_r.Draw();
                    if abs(V2[2]-0.5)<abs(V2[3]-0.5): self.exponent = V2[2]; 
                    if abs(V2[2]-0.5)>abs(V2[3]-0.5): self.exponent = V2[3];
                #======= GET MEANS and COVARIANCE ============
                #Get_MuVec_CovMat
                #Get_MeanVec_CovMat(h_Signal1 ,h_TotalMC);
                ###V2 = OptimalCut(h_TotalMC, h_Signal2, NORM_s1);
                ###sig2_Leftcut =TGaxis( V2[2], MinY, V2[2],MaxY, 0,0,0, "L"); sig2_Leftcut.SetLineColor(kGreen+2);  sig2_Leftcut.SetLineWidth(3);  sig2_Leftcut.Draw();
                ###sig2_Rightcut=TGaxis( V2[3], MinY, V2[3],MaxY, 0,0,0, "L"); sig2_Rightcut.SetLineColor(kGreen+2); sig2_Rightcut.SetLineWidth(3); sig2_Rightcut.Draw();
                ###PrintValus2 = TLatex( h_TotalMC.GetBinLowEdge(2), MaxY*0.5,"    r "+V2[4]+"%   e "+V2[5]+"%   s "+V2[6]+"%"); PrintValus2.SetTextColor(kGreen+2); PrintValus2.SetTextSize(0.13); PrintValus2.Draw();
                ###ar_l = TArrow( V2[2] ,MaxY/2, V2[2]+length, MaxY/2, 0.03 , "|>"); ar_l.SetLineWidth(3); ar_l.SetFillColor(kGreen+2); ar_l.SetLineColor(kGreen+2); ar_l.SetAngle(30); ar_l.Draw();
                ###ar_r = TArrow( V2[3] ,MaxY/2, V2[3]-length, MaxY/2, 0.03 , "|>"); ar_r.SetLineWidth(3); ar_r.SetFillColor(kGreen+2); ar_r.SetLineColor(kGreen+2); ar_r.SetAngle(30); ar_r.Draw();
                ###if abs(V2[2]-0.5)<abs(V2[3]-0.5): self.exponent = V2[2];
                ###if abs(V2[2]-0.5)>abs(V2[3]-0.5): self.exponent = V2[3];

            fPads2.RedrawAxis(); fPads2.Update();
            fPads1.RedrawAxis(); fPads1.Update();

        if MODE in ["COMP"]:
            print " ---> bottom pad"; 
            fPads2.cd();
            h_Ratio=TH1D("h_Ratio",""+";%s;%s"%(xtitle,ytitle),nbin,min,max)
            h_Ratio.SetLineColor(1)
            h_Ratio.SetLineWidth(2)
            h_Ratio.SetMarkerColor(1)
            h_Ratio.SetMarkerStyle(8)
            h_Ratio.SetMarkerSize(0.7)
            h_Ratio.SetStats(0)
            maxY=2
            axis1=TGaxis( min,1,max,1, 0,0,0, "L"); axis1.SetLineColor(1); axis1.SetLineWidth(1);  #axis1->SetLabelColor(16);
            for i in range(1,h_Ratio.GetNbinsX()+1,1):
                D  = h_DMR.GetBinContent(i);    eD = h_DMR.GetBinError(i);
                if D==0: eD=0.92;
                B  = h_QCD.GetBinContent(i); eB = h_QCD.GetBinError(i);
                #if options.COMP[2]=="t":
                #      B  = h_T.GetBinContent(i); eB = h_T.GetBinError(i);
                if B<0.1 and eB>=B : eB=0.92; Err= 0.;
                if B!=0.        :Err=TMath.Sqrt( (eD*eD)/(B*B)  +(D*D*eB*eB)/(B*B*B*B)     ); h_Ratio.SetBinContent(i, D/B   );  h_Ratio.SetBinError(i, Err); #print i,")",h_Ratio.GetNbinsX()+1,")   data:",D," pm ",eD,"     Bkg:",B," pm ",eB,"   R:",D/B," pm ", Err
                if B==0.        :Err=TMath.Sqrt( (eD*eD)/(eB*eB)+(D*D*eB*eB)/(eB*eB*eB*eB) ); h_Ratio.SetBinContent(i, D/0.92);  h_Ratio.SetBinError(i, Err);
                if D==0 and B==0:                                                             h_Ratio.SetBinContent(i, -1);      h_Ratio.SetBinError(i, 0  );
                if h_Ratio.GetBinContent(i)>maxY:h_Ratio.SetBinContent(i, maxY); ### To visualise the points above axis... #h_Ratio.Fit("pol1");
            h_Ratio2=TH1D("h_Ratio2",""+";%s;%s"%(xtitle,ytitle),nbin,min,max)
            h_Ratio2.SetLineColor(2)
            h_Ratio2.SetLineWidth(3)
            h_Ratio2.SetMarkerColor(2)
            h_Ratio2.SetMarkerStyle(8);
            h_Ratio2.SetMarkerSize(0);
            h_Ratio2.GetYaxis().SetRangeUser( 0 , maxY )
            h_Ratio2.GetYaxis().SetNdivisions(504,0)
            h_Ratio2.GetYaxis().SetTitle("Ratios   ")
            h_Ratio2.GetYaxis().SetTitleOffset(0.35)
            h_Ratio2.GetYaxis().SetTitleSize(0.13)
            h_Ratio2.GetYaxis().SetTitleSize(0.13)
            h_Ratio2.GetYaxis().SetLabelSize(0.11)
            h_Ratio2.GetXaxis().SetLabelSize(0.1)
            h_Ratio2.GetXaxis().SetTitleOffset(1.0)
            h_Ratio2.GetXaxis().SetTitleSize(0.1)
            h_Ratio2.SetStats(0)

            for i in range(1,h_Ratio2.GetNbinsX()+1,1):
                D  = h_QCD_S.GetBinContent(i);    eD = h_QCD_S.GetBinError(i);
                B  = h_QCD.GetBinContent(i); eB = h_QCD.GetBinError(i);
            #if options.COMP[2]=="t":
            #    D  = h_T_S.GetBinContent(i);    eD = h_T_S.GetBinError(i);
            #    B  = h_T.GetBinContent(i); eB = h_T.GetBinError(i);
                if D==0: eD=0.92;
                if B<0.1 and eB>=B : eB=0.92; Err= 0.;
                if B!=0.        :Err=TMath.Sqrt( (eD*eD)/(B*B)  +(D*D*eB*eB)/(B*B*B*B)     ); h_Ratio2.SetBinContent(i, D/B   );  h_Ratio2.SetBinError(i, Err);
                if B==0.        :Err=TMath.Sqrt( (eD*eD)/(eB*eB)+(D*D*eB*eB)/(eB*eB*eB*eB) ); h_Ratio2.SetBinContent(i, D/0.92);  h_Ratio2.SetBinError(i, Err);
                if D==0 and B==0:                                                             h_Ratio2.SetBinContent(i, -1);      h_Ratio2.SetBinError(i, 0  );
                if h_Ratio2.GetBinContent(i)>maxY:h_Ratio2.SetBinContent(i, maxY); ### To visualise the points above axis... #h_Ratio.Fit("pol1");
            axis1=TGaxis( min,1,max,1, 0,0,0, "L")
            axis1.SetLineColor(1)
            axis1.SetLineWidth(1)
            h_Ratio2.Draw("hist e0")
            h_Ratio.Draw("same e0")
            axis1.Draw()
            theLeg2 = TLegend(0.23, 0.23, 0.4, 0.5,"","NDC"); theLeg2.SetName("theLegend2"); theLeg2.SetBorderSize(0); theLeg2.SetLineColor(0); theLeg2.SetFillColor(0);theLeg2.SetFillStyle(0); theLeg2.SetLineWidth(0); theLeg2.SetLineStyle(0); theLeg2.SetTextFont(42);theLeg2.SetTextSize(.1);
            theLeg2.SetFillColor(0);theLeg2.SetBorderSize(0);theLeg2.SetLineColor(0);theLeg2.SetLineWidth(0);theLeg2.SetLineStyle(0);theLeg2.SetTextFont(42);
            theLeg2.AddEntry(h_Ratio, "[Data-Rest]_CR / QCD_CR","ape");
            theLeg2.AddEntry(h_Ratio2, "QCD_SR / QCD_CR","ape");
            theLeg2.Draw();


        #============= THE LEGEND SESSION =======================
        if MODE in ["MC","MCvsDATA"]:
            theLeg = TLegend(0.48, 0.55, 0.9, 0.9, "", "NDC");theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(.05);
            theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42);#theLeg.SetNColumns(2);
            if MODE=="MCvsDATA":theLeg.AddEntry(h_data, "Data "+options.y,"ep");
            theLeg.AddEntry(h_QCD    , "QCD"             ,"F");
            theLeg.AddEntry(h_WJets  , "W+jets "         ,"F");
            theLeg.AddEntry(h_TTbar  , "t#bar{t}"        ,"F");
            theLeg.AddEntry(h_STop   , "single t"        ,"F");
            theLeg.AddEntry(h_Rest   , "Rest (VV,Z+jets)","F"); #end1=""; end2="";
            if REGION[:3] in ["PS2","CR1","CR2","CR3","SR1","SR2","SR3"]:  theLeg.AddEntry(h_Signal1,"(2000,200) #times %s "%(NORM_s1),"L");
            if REGION[:3] in ["PS3","CR4","SR4"]:                          theLeg.AddEntry(h_Signal2,"(3000,1500) #times %s "%(NORM_s2),"L");
            #theLeg.AddEntry(h_Signal2,"(3000, 1500) #times %s "%(NORM_s2),"L");
            theLeg.SetY1NDC(0.9-0.08*6-0.005);
            theLeg.SetY1(theLeg.GetY1NDC()); fPads1.cd(); theLeg.Draw(); #theLeg.AddEntry(gr_MCStat, "Sys.","F");
            #============ SET MAX Y-AXIS FOR PLOTS ==================
            histsigmax = TMath.Max( h_Signal1.GetMaximum(), h_Signal2.GetMaximum() );            histsigmin = TMath.Min( h_Signal1.GetMinimum(), h_Signal2.GetMinimum() );
            #histsigmax = TMath.Max( histsigmax, h_Signal3.GetMaximum() );
            #histsigmax = TMath.Max( histsigmax, h_Signal4.GetMaximum() );
            histsigmax = TMath.Max( histsigmax, h_data.GetMaximum() );            histsigmin = TMath.Min( histsigmin, h_data.GetMinimum() );
            histsigmax = TMath.Max( histsigmax, h_TotalMC.GetMaximum() );         histsigmin = TMath.Min( histsigmin, h_TotalMC.GetMinimum() );
            h_Signal1.GetYaxis().SetRangeUser(0, histsigmax*1.3 );
            h_TotalMC.GetYaxis().SetRangeUser(0, histsigmax*1.3 );
            h_data.GetYaxis().SetRangeUser(   0, histsigmax*1.3 );
            if logy == 1:
                if histsigmin<=0: histsigmin=1;
                h_Signal1.GetYaxis().SetRangeUser(0.5*histsigmin, histsigmax*100. );
                h_TotalMC.GetYaxis().SetRangeUser(0.5*histsigmin, histsigmax*100. );
                h_data.GetYaxis().SetRangeUser(   0.5*histsigmin, histsigmax*100. );
        if MODE in ["DECO"]:
            theLeg = TLegend(0.5, 0.75, 0.75, 0.92, "", "NDC");theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(.04);
            theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42);#theLeg.SetNColumns(2);
            theLeg.AddEntry(h_TotalMC,"Total BKG shape","L");
            theLeg.AddEntry(h_Signal1,"Total signal (2000,200) GeV","L");
            #theLeg.AddEntry(h_TotalS1,"Signal as sum of comp.","L");
            print "Non-Closure of the sum of components:",(h_TotalS1.Integral()-h_Signal1.Integral())/(h_Signal1.Integral()+0.000001);
            #theLeg.AddEntry(h_Signal1_R4q ,"R^{4q}  %.0f"%(  h_Signal1_R4q.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            #theLeg.AddEntry(h_Signal1_R3q ,"R^{3q}  %.0f"%(  h_Signal1_R3q.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            #theLeg.AddEntry(h_Signal1_R2q ,"R^{2q}  %.0f"%(  h_Signal1_R2q.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            #theLeg.AddEntry(h_Signal1_W   ,"W       %.0f"%(    h_Signal1_W.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            #theLeg.AddEntry(h_Signal1_Rlqq,"R^{lqq} %.0f"%( h_Signal1_Rlqq.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            #theLeg.AddEntry(h_Signal1_Rlq ,"R^{lq}  %.0f"%(  h_Signal1_Rlq.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            #theLeg.AddEntry(h_Signal1_g   ,"g       %.0f"%(    h_Signal1_g.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            #theLeg.AddEntry(h_Signal1_Rest,"Rest    %.0f"%( h_Signal1_Rest.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            theLeg.AddEntry(h_Signal1_Neut,"Neutrino component %.0f"%(h_Signal1_Neut.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            #theLeg.AddEntry(h_Signal1_qlea,"q leakage %.0f"%(h_Signal1_qlea.Integral()/(h_TotalS1.Integral()+0.00001)*100)+"%" ,"L");
            theLeg.SetY1NDC(0.9-0.08*6-0.005);
            theLeg.SetY1(theLeg.GetY1NDC()); theLeg.Draw();

        if MODE in ["COMP"]:
            print " ---> test"; 
            theLeg = TLegend(0.45, 0.8, 0.6, 0.6, "", "NDC");theLeg.SetName("theLegend"); theLeg.SetBorderSize(0); theLeg.SetLineColor(0); theLeg.SetFillColor(0);theLeg.SetFillStyle(0); theLeg.SetLineWidth(0); theLeg.SetLineStyle(0); theLeg.SetTextFont(42);theLeg.SetTextSize(.05);
            theLeg.SetFillColor(0);theLeg.SetBorderSize(0);theLeg.SetLineColor(0);theLeg.SetLineWidth(0);theLeg.SetLineStyle(0);theLeg.SetTextFont(42);#theLeg.SetNColumns(2);
            theLeg.AddEntry(h_QCD_S, "QCD at S"+options.REGION[1:]+" scaled by "+str(round(NORM,2))  ,"L");
            theLeg.AddEntry(h_QCD  , "QCD at C"+options.REGION[1:3]+" norm. at S"+options.REGION[1:] ,"L");
            theLeg.AddEntry(h_DMR  , "Data-Rest_MC norm. at S"+options.REGION[1:],"L");
            histsigmax = TMath.Max(h_QCD_S.GetMaximum(), h_QCD.GetMaximum() );
            histsigmax = TMath.Max( histsigmax, h_DMR.GetMaximum() );
            h_QCD_S.GetYaxis().SetRangeUser(0, histsigmax*1.3 );
            h_QCD.GetYaxis().SetRangeUser(  0, histsigmax*1.3 );
            h_DMR.GetYaxis().SetRangeUser(  0, histsigmax*1.3 );
            theLeg.SetY1NDC(0.9-0.08*6-0.005);
            theLeg.SetY1(theLeg.GetY1NDC()); fPads1.cd(); theLeg.Draw(); #theLeg.AddEntry(gr_MCStat, "Sys.","F");

            
        #============ SAVE PLOTS IN A DIRECTORY ============================
        extension   = "";
        if tag    !=  "": extension = extension + "_"+tag;
        if logy         : extension = extension + "_log";
        if options.FBT  : extension = extension + "_FBT";
        if SFs          : extension = extension + "_SFsCorr";
        #----------- Rename variables to a shorter name -----------------
        for c in [".","/","(",")","[","]","*","+",">","<"," ","=",",",":","deep","dnn","Decorr","jetAK8puppi","ass_tag","t_tag","_tag","TMath","Cos","Sin"]:variable=variable.replace(c,"_");
        for c in ["__","___","____","_____","______","_"]:variable=variable.replace(c,"");
        #for c in ["__","___","____","_____","______","_"]:variable=variable.replace(c,"_");
        #----------------- Save and open the plot -----------------------
        #Name=options.TAG+REGION+"_"+MODE+"_"+variable+"_"+options.y+extension+".png";
        #Name=tag+REGION+"_"+MODE+"__"+options.y+extension+".png"
        Name=variable+"_"+tag+REGION+"_"+MODE+"_"+options.y+extension+".png"
        #Name=V1[6]+"_"+options.TAG+REGION+"_"+MODE+"_"+variable+"_"+options.y+extension+".png"        

        file=TString(Name); 
        canvas_controlplot.SaveAs( file.Data() )
        os.system("display %s &"%(Name) ); print " --> display %s &"%(Name);


        if options.piechart:#============== PIE CHARTS =================================
            num_events, colors = array('d'), array('i');   gStyle.SetOptStat(000000000);            
            piecanvas=TCanvas("PIES","PIES",400,400);  piecanvas.SetTickx(1); piecanvas.SetTicky(1); piecanvas.SetRightMargin(-0.5); piecanvas.SetTopMargin(-0.5); piecanvas.SetLeftMargin(-0.5); piecanvas.SetBottomMargin(-0.5);gStyle.SetOptStat(000000000);
            num_events.append(h_QCD.Integral());   num_events.append(h_WJets.Integral());  num_events.append(h_TTbar.Integral()); num_events.append(h_STop.Integral()); num_events.append(h_Rest.Integral());
            colors.append(self.color_palet["QCD"]);   colors.append(self.color_palet["WJets"]);  colors.append(self.color_palet["TTbar"]); colors.append(self.color_palet["STop"]); colors.append(self.color_palet["Rest"]);
            pieplot=TPie("PIE","",5,num_events,colors);
            pieplot.SetEntryLabel(0,"QCD");pieplot.SetEntryLabel(1,"W+jets"); pieplot.SetEntryLabel(2,"t#bar{t}");pieplot.SetEntryLabel(3,"single t");pieplot.SetEntryLabel(4,"Rest VV,Z+jets");
            pieplot.SetTextSize(.045); pieplot.SetAngularOffset(30); pieplot.SetLabelFormat("%val (%perc) %txt"); pieplot.SetRadius(.4); pieplot.SetLabelsOffset(-.33); piecanvas.cd(1);pieplot.Draw("nol <");RegionTxt.Draw();
            file2=TString("PIE_%s_%s.png"%(REGION,options.y)); piecanvas.SaveAs(file2.Data());
            os.system("display PIE_%s_%s.png &"%(REGION,options.y) );




################# Main Code ################################
def Draw_Control_Plot( channel ) :
    if channel in ["had"]           : Instance_ANALYSIS = ANALYSIS( channel ); Instance_ANALYSIS.DefineSelection_0lep();
    #if channel in ["mu","el","lep"] : Instance_ANALYSIS = ANALYSIS( channel ); Instance_ANALYSIS.DefineSelection_1lep();

if __name__ == '__main__' : 
    Beginning = strftime("%H:%M:%S",gmtime())
    print '\n----RUN--------------channel:[',options.channel,']----------Region:[',options.REGION,']----------------[',Beginning,']--------'
    Draw_Control_Plot( options.channel );
    Finishing = strftime("%H:%M:%S",gmtime());
    #========== CALCULATE DURATION OF THE RUN ===========
    MIN=int(Finishing[3:5])-int(Beginning[3:5]); SEC=int(Finishing[6:8])-int(Beginning[6:8]); 
    if SEC<0 and MIN>0 : SEC=60+SEC; MIN=MIN-1;
    if SEC>0 and MIN<0 : MIN=60+MIN;
    if SEC<0 and MIN<0 : SEC=60+SEC; MIN=60+MIN-1;
    print '----END-----------------------------------------------[time:',Finishing,', duration:',MIN,':',SEC,']---------\n'#!usr/bin/env python



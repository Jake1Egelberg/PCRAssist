
#Package install

package.list<-c("tcltk",
                "stringr")
for(i in 1:length(package.list)){
  tryCatch(find.package(package.list[i]),
           error = function(e) install.packages(package.list[i],repos="http://lib.stat.cmu.edu/R/CRAN/"))
}

library(tcltk)
library(stringr)

#Load variables

primers_in_master_var<-tclVar(0)
enzyme_conc_X_var<-tclVar(2)
rxt_volume_ul_var<-tclVar(25)
rxt_num_var<-tclVar(3)

init_f_r_prm_concs<-tclVar("10,10")
fin_f_r_prm_concs<-tclVar("0.4,0.4")

sample_volume_uL_var<-tclVar(100)
dH2O_volume<-NA

sel_dna_var<-tclVar("Genomic")

#Create font
large<-tkfont.create(size=15)

#Display output function
create_display<-function(){
  display<-tktoplevel()
  tkwm.title(display,"Protocol")
  tkwm.geometry(display,"620x280+350+100")
  frame<-tkframe(display)
  tkgrid(frame,sticky="w",padx=10,pady=10)
  prot<-tklabel(frame,text=protocol,justify="left",font=large)
  tkgrid(prot)
  tkwait.window(display)
}

#Load variables
load_vars<-function(){
  .GlobalEnv$primers_in_master_val<-tclvalue(primers_in_master_var)
  if(primers_in_master_val=="1"){
    .GlobalEnv$primers_in_master<-TRUE
  } else{
    .GlobalEnv$primers_in_master<-FALSE
  }
  
  .GlobalEnv$enzyme_conc_X<-as.numeric(tclvalue(enzyme_conc_X_var))
  .GlobalEnv$rxt_volume_ul<-as.numeric(tclvalue(rxt_volume_ul_var))
  .GlobalEnv$rxt_num<-as.numeric(tclvalue(rxt_num_var))
  
  .GlobalEnv$init_f_r_prm_concs_vals<-tclvalue(init_f_r_prm_concs)
  init_splits<-str_split(init_f_r_prm_concs_vals,",")[[1]]
  .GlobalEnv$f_primer_conc_uM<-as.numeric(init_splits[1])
  .GlobalEnv$r_primer_conc_uM<-as.numeric(init_splits[2])
  
  .GlobalEnv$fin_f_r_prm_concs_vals<-tclvalue(fin_f_r_prm_concs)
  fin_splits<-str_split(fin_f_r_prm_concs_vals,",")[[1]]
  .GlobalEnv$desired_f_conc_uM<-as.numeric(fin_splits[1])
  .GlobalEnv$desired_r_conc_uM<-as.numeric(fin_splits[2])
  
  .GlobalEnv$sample_volume_uL<-as.numeric(tclvalue(sample_volume_uL_var))
  .GlobalEnv$dH2O_volume<-NA
  
  .GlobalEnv$sel_dna_val<-as.character(tclvalue(sel_dna_var))
  
}

#Do math
calculate_metrics<-function(){
  
  load_vars()
  
  if(sel_dna_val=="Genomic"){
    #Add 35 ng DNA to master mix
    .GlobalEnv$amt_added<-35
    .GlobalEnv$sample_volume_to_add<-(amt_added/sample_volume_uL)
  } else if(sel_dna_val=="Vector"){
    #Add 1 ng DNA to master mix
    .GlobalEnv$amt_added<-1
    .GlobalEnv$sample_volume_to_add<-(amt_added/sample_volume_uL)
  } else if(sel_dna_val=="Colony"){
    .GlobalEnv$amt_added<-"?"
    .GlobalEnv$sample_volume_to_add<-2
  }
  
  if(primers_in_master==FALSE){
    
    #GoTaq is 2X concentrated
    gotaq_volume_2x<-rxt_volume_ul/enzyme_conc_X
    
    #Dilute primers directly into PCR tubes from working stock
    f_primer_volume<-(desired_f_conc_uM/f_primer_conc_uM)*rxt_volume_ul
    r_primer_volume<-(desired_r_conc_uM/r_primer_conc_uM)*rxt_volume_ul
    
    #Get remaining dH2O per tube
    dH2O_volume<-rxt_volume_ul-(gotaq_volume_2x+f_primer_volume+r_primer_volume+sample_volume_to_add)
    
    #Get GoTAQ & dH2O per tube
    non_prim_vol<-dH2O_volume+gotaq_volume_2x
    
    #Calculate GoTAQ and dH2O in master mix
    gotaq_master<-gotaq_volume_2x*(rxt_num+1)
    dH2O_master<-dH2O_volume*(rxt_num+1)
    
    .GlobalEnv$protocol<-paste("PRIMERS INDIVIDUAL\n\n1) Mix ",round(gotaq_master,digits=1)," uL ",round(enzyme_conc_X,digits=1),"X polymerase with ",round(dH2O_master,digits=1)," uL dH2O for master mix. \n\n2) Aliquot ",round(non_prim_vol,digits=1)," uL of master mix into each tube. \n\n3) Add ",round(f_primer_volume,digits=1)," uL of forward primer stock to each tube. \n\n4) Add ",round(r_primer_volume,digits=1)," uL of reverse primer stock to each tube. \n\n5) Add ",round(sample_volume_to_add,digits=2)," uL (",amt_added," ng) of ",paste(sel_dna_val," DNA",sep="")," to each tube.",sep="")
    
    create_display()
    
  } else if(primers_in_master==TRUE){
    
    #Master mix volume per tube
    #GoTaq, dH2O, primers
    mm_vol<-rxt_volume_ul-sample_volume_to_add
    
    #Total master mix to make
    total_mm<-mm_vol*(rxt_num+1)
    
    #Primer concentrations in master mix
    f_primer_conc<-(rxt_volume_ul/mm_vol)*desired_f_conc_uM
    r_primer_conc<-(rxt_volume_ul/mm_vol)*desired_r_conc_uM
    #Primer volumes to add
    f_primer_volume<-(f_primer_conc/f_primer_conc_uM)*total_mm
    r_primer_volume<-(r_primer_conc/r_primer_conc_uM)*total_mm
    
    #GoTaq conc in master mix
    gotaq_conc<-(rxt_volume_ul/mm_vol)*1
    #GoTaq volume in master mix
    gotaq_vol<-(gotaq_conc/enzyme_conc_X)*total_mm
    
    #dH2O in mm
    dH2O_mm<-total_mm-(gotaq_vol+f_primer_volume+r_primer_volume)
    
    .GlobalEnv$protocol<-paste("PRIMERS IN MASTER MIX\n\n1) Mix ",round(gotaq_vol,digits=1)," uL ",round(enzyme_conc_X,digits=1),"X polymerase with ",round(dH2O_mm,digits=1)," uL dH2O for master mix. \n\n2) Add ",round(f_primer_volume,digits=1)," uL of forward primer to master mix. \n\n3) Add ",round(r_primer_volume,digits=1)," uL of reverse primer stock to master mix. \n\n4) Aliquot ",round(mm_vol,digits=1)," uL of master mix into each tube. \n\n5) Add ",round(sample_volume_to_add,digits=2)," uL (",amt_added," ng) of ",paste(sel_dna_val," DNA",sep="")," to each tube.",sep="")
    
    create_display()
    
  }
  
}

#Create gui
gui<-tktoplevel()
tkwm.title(gui,"PCR Assist")
tkwm.geometry(gui,"240x430+100+100")

frm<-tkframe(gui)
tkgrid(frm,padx=20,column=1,row=1,columnspan=2)
ttl<-tklabel(frm,text="PCR Assist",font=tkfont.create(size=15,weight="bold"))
tkgrid(ttl,column=1,columnspan=2,pady=5)

rxt_parms<-tklabel(frm,text="Reaction Parameters",font=tkfont.create(weight="bold",size=10))
tkgrid(rxt_parms,row=1,column=1,columnspan=2)

rx_vol_lbl<-tklabel(frm,text="Reaction volume (uL)")
tkgrid(rx_vol_lbl,pady=5,sticky="w",row=2,column=1)
rx_vol<-tkentry(frm,textvariable=rxt_volume_ul_var,justify="center",width=11)
tkgrid(rx_vol,sticky="e",pady=5,row=2,column=2)
tkbind(rx_vol,"<Return>",calculate_metrics)
rx_num_lbl<-tklabel(frm,text="Reaction #")
tkgrid(rx_num_lbl,pady=5,sticky="w",row=3,column=1)
rx_num<-tkentry(frm,textvariable=rxt_num_var,justify="center",width=11)
tkgrid(rx_num,sticky="e",pady=5,row=3,column=2)
tkbind(rx_num,"<Return>",calculate_metrics)
polym_conc_lbl<-tklabel(frm,text="[Polymerase] (X)")
tkgrid(polym_conc_lbl,pady=5,sticky="w",row=4,column=1)
polym_conc<-tkentry(frm,textvariable=enzyme_conc_X_var,justify="center",width=11)
tkgrid(polym_conc,sticky="e",pady=5,row=4,column=2)
tkbind(polym_conc,"<Return>",calculate_metrics)

samp_parms<-tklabel(frm,text="Sample Parameters",font=tkfont.create(weight="bold",size=10))
tkgrid(samp_parms,pady=5,row=5,column=1,columnspan=2)

samp_vol<-tklabel(frm,text="[DNA Sample] (ng/uL)")
tkgrid(samp_vol,pady=5,sticky="w",row=6,column=1)
samp_vol_entr<-tkentry(frm,textvariable=sample_volume_uL_var,justify="center",width=11)
tkgrid(samp_vol_entr,sticky="e",pady=5,row=6,column=2)
tkbind(samp_vol_entr,"<Return>",calculate_metrics)
dna_type<-tklabel(frm,text="Sample type")
tkgrid(dna_type,pady=5,sticky="w",row=7,column=1)
sel_dna<-ttkcombobox(frm,values=c("Genomic","Vector","Colony"),width=8,justify="center",textvariable=sel_dna_var)
tkgrid(sel_dna,sticky="w",pady=5,row=7,column=2)

prim_parms<-tklabel(frm,text="Primer Parameters",font=tkfont.create(weight="bold",size=10))
tkgrid(prim_parms,pady=5,row=8,column=1,columnspan=2)

primers_in_master_but<-tkcheckbutton(frm,text="Primers in master mix?",variable=primers_in_master_var)
tkgrid(primers_in_master_but,pady=2,row=9,column=1,columnspan=2)
init_primer_concs<-tklabel(frm,text="Initial F,R [Primer] (uM)")
tkgrid(init_primer_concs,pady=5,sticky="w",row=10,column=1)
init_primer_concs_entr<-tkentry(frm,textvariable=init_f_r_prm_concs,justify="center",width=11)
tkgrid(init_primer_concs_entr,sticky="e",row=10,pady=5,column=2)
tkbind(init_primer_concs_entr,"<Return>",calculate_metrics)
des_primer_concs<-tklabel(frm,text="Desired F,R [Primer] (uM)")
tkgrid(des_primer_concs,pady=5,sticky="w",row=11,column=1)
des_primer_concs_entr<-tkentry(frm,textvariable=fin_f_r_prm_concs,justify="center",width=11)
tkgrid(des_primer_concs_entr,sticky="e",row=11,pady=5,column=2)
tkbind(des_primer_concs_entr,"<Return>",calculate_metrics)

gen_prot<-tkbutton(frm,text="Generate protocol",command=calculate_metrics)
tkgrid(gen_prot,pady=10,row=12,column=1,columnspan=2)

tkwait.window(gui)


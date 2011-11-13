'''

Implementation of mu-mu-tau channel analysis

'''

import ROOT
import os
import sys
import logging
import FinalStateAnalysis.PatTools.data as data_tool
from FinalStateAnalysis.Utilities.AnalysisPlotter import styling,samplestyles

logging.basicConfig(filename='mmtAnalysis.log',level=logging.DEBUG, filemode='w')
log = logging.getLogger("mmtChannel")

ROOT.gROOT.SetBatch(True)

# Define the fake rate versus muon pt
FAKE_RATE_0 = 0.00435
FAKE_RATE_1 = 4.247
FAKE_RATE_2 = -0.1967

FAKE_RATE = '(%0.4f + %0.4f*TMath::Exp(%0.4f * (VAR) ) )' % (
    FAKE_RATE_0, FAKE_RATE_1, FAKE_RATE_2)

FR_X = 'Muon2Pt + Muon2Pt*Muon2_MuRelIso'
FAKE_RATE = FAKE_RATE.replace('VAR', FR_X)
#FR_WEIGHT = '((%s)/(1-%s))' % (FAKE_RATE, FAKE_RATE)
FR_WEIGHT = '(%s)' % FAKE_RATE

base_selection = [
    #'(run < 5 && Muon1Pt > 13.5 || Muon1Pt > 13.365)',
    #'(run < 5 && Muon2Pt > 9 || Muon2Pt > 8.91)',
    'Muon1Pt > 18',
    'Muon2Pt > 9',
    'Muon1AbsEta < 2.1',
    'Muon2AbsEta < 2.1',
    'Muon1_MuRelIso < 0.3',
    'Muon1_MuID_WWID > 0.5',
    #'NIsoMuonsPt5_Nmuons < 0.5',
    #'DoubleMus_HLT > 0.5 ',
    #'NIsoTausPt20_NIsoTaus < 0.5',
    #'TauPt > 20',
    'Tau_LooseHPS > 0.5',
    'Muon1Charge*Muon2Charge > 0',
    'Muon2_NPixHits > 0',
]

passes_ht = [
    'VisFinalState_Ht > 80'
]

passes_vtx = [
    'vtxChi2/vtxNDOF < 15'
]

bkg_enriched = [
    '(Muon2_MuID_WWID < 0.5 || Muon2_MuRelIso > 0.3)'
]

final_selection = [
    'Muon2_MuID_WWID > 0.5',
    'Muon2_MuRelIso < 0.3'
]

variables = {
#    'DiMuonMass' : ('Muon1_Muon2_Mass', 'M_{#mu#mu}', [100, 0, 300],),
    'MuTauMass' : ('Muon2_Tau_Mass', 'M_{#mu#tau}', [60, 0, 300],),
    'Muon1_MtToMET' : ('Muon1_MtToMET', 'M_{T} #mu(1)-#tau', [60, 0, 300],),
#    'Muon2_MtToMET' : ('Muon2_MtToMET', 'M_{T} #mu(2)-#tau', [100, 0, 300],),
    'vtxChi2NODF' : ('vtxChi2/vtxNDOF', 'Vertex #chi^{2}/NODF', [100, 0, 30],),
#    'MET' : ('METPt', 'MET', [100, 0, 200]),
    'HT' : ('VisFinalState_Ht', 'H_{T}', [60, 0, 300]),
}

selections = {
    'base_line' : {
        'select' : base_selection,
        'title' : "Base Selection",
        'vars' : [
            'DiMuonMass',
            'MuTauMass',
            'Muon1_MtToMET',
            'MET',
            'HT',
            'vtxChi2NODF',
        ],
    },
    'with_ht' : {
        'select' : base_selection + passes_ht,
        'title' : "Base Selection + H_{T}",
        'vars' : [
            'DiMuonMass',
            'MuTauMass',
            'Muon1_MtToMET',
            'vtxChi2NODF',
        ],
    },
    'final' : {
        'title' : "Final Selection",
        'select' : base_selection + passes_ht + passes_vtx,
        'vars' : [
            'DiMuonMass',
            'MuTauMass',
            'Muon1_MtToMET',
            'vtxChi2NODF',
        ],
    },
}

# Samples we aren't interested in
skips = ['MuEG', 'DoubleEl', 'EM', 'v1_d', ]
int_lumi = 3900
samples, plotter = data_tool.build_data(
    'VH', '2011-11-07-v1-WHAnalyze', 'scratch_results',
    int_lumi, skips, count='emt/skimCounter')

canvas = ROOT.TCanvas("basdf", "aasdf", 800, 600)
def saveplot(filename):
    # Save the current canvas
    filetype = '.pdf'
    #legend.Draw()
    canvas.SetLogy(False)
    canvas.Update()
    canvas.Print(os.path.join(
        "plots", 'mmtChannel', filename + filetype))
    canvas.SetLogy(True)
    canvas.Update()
    canvas.Print(os.path.join(
        "plots", 'mmtChannel', filename + '_log' + filetype))

################################################################################
###  Control plot: check integrated lumi   #####################################
################################################################################

# Version with weights applied
lumi = plotter.get_histogram( 'data_DoubleMu', '/mmt/intLumi',)
lumi.Draw()
saveplot('intlumi')

for selection, selection_info in selections.iteritems():
    log.info("Doing " + selection)

    log.info("Plotting distribution of FR weights")

    fr_selection = selection_info['select'] + bkg_enriched
    the_final_selection = selection_info['select'] + final_selection

    # Version with weights applied
    plotter.register_tree(
        selection + '_fr_weights',
        '/mmt/final/Ntuple',
        FR_WEIGHT,
        ' && '.join(fr_selection),
        #w = '(%s)' % FR_WEIGHT,
        binning = [200, 0, 1],
        include = ['*data*'],
    )

    weight_histo = plotter.get_histogram(
        'data_DoubleMu',
        '/mmt/final/Ntuple:' + selection + '_fr_weights',
        show_overflows = True,
    )
    weight_histo.Draw()
    saveplot(selection + '_fr_weights')

    for var in selection_info['vars']:
        if var not in variables:
            print "Skipping", var
            continue
        draw_str, x_title, binning = variables[var]

        plotter.register_tree(
            selection + var + '_bkg',
            '/mmt/final/Ntuple',
            draw_str,
            ' && '.join(fr_selection),
            w = '(puWeight_3bx_S42011AB178078)',
            binning = binning,
            include = ['*'],
        )

        # Version with weights applied
        plotter.register_tree(
            selection + var + '_bkg_fr',
            '/mmt/final/Ntuple',
            draw_str,
            ' && '.join(fr_selection),
            w = '(puWeight_3bx_S42011AB178078)*(%s)' % FR_WEIGHT,
            binning = binning,
            include = ['*'],
        )

        plotter.register_tree(
            selection + var + '_fin',
            '/mmt/final/Ntuple',
            draw_str,
            ' && '.join(the_final_selection),
            w = 'puWeight_3bx_S42011AB178078',
            binning = binning,
            include = ['*'],
        )

        ##########################################
        # Compute results using MC in final region
        ##########################################
        rebin = 5

        histo_name = selection + var + '_bkg'

        stack = plotter.build_stack(
            '/mmt/final/Ntuple:' + histo_name,
            include = ['*'],
            exclude = ['data*'],
            rebin = rebin, show_overflows=True,
        )

        data = plotter.get_histogram(
            'data_DoubleMu',
            '/mmt/final/Ntuple:' + histo_name,
            rebin = rebin, show_overflows=True,
        )

        stack.Draw()
        data.Draw('same,pe')
        saveplot(histo_name)

        histo_name = selection + var + '_fin'

        stack = plotter.build_stack(
            '/mmt/final/Ntuple:' + histo_name,
            include = ['*'],
            exclude = ['data*'],
            rebin = rebin, show_overflows=True,
        )

        data = plotter.get_histogram(
            'data_DoubleMu',
            '/mmt/final/Ntuple:' + histo_name,
            rebin = rebin, show_overflows=True,
        )

        stack.Draw()
        data.Draw('same,pe')
        stack.SetMaximum(max(stack.GetMaximum(), data.GetMaximum())*1.5)
        stack.GetXaxis().SetTitle(x_title)
        saveplot(histo_name)

        ##########################################
        # Compute results using FR metho
        ##########################################

        # The final selected events
        data = plotter.get_histogram(
            'data_DoubleMu',
            '/mmt/final/Ntuple:' + selection + var + '_fin',
            rebin = rebin, show_overflows = True
        )

        # The background prediction
        data_bkg_fr = plotter.get_histogram(
            'data_DoubleMu',
            '/mmt/final/Ntuple:' + selection + var + '_bkg_fr',
            rebin = rebin, show_overflows = True
        )

        corrected_mc = ['ZZ', 'WZ']
        corrected_mc_histos = []

        # Now we need to correct the backgrounds which aren't correctly
        # estimated by the FR method.
        for to_correct in corrected_mc:
            mc_final = plotter.get_histogram(
                to_correct,
                '/mmt/final/Ntuple:' + selection + var + '_fin',
                rebin = rebin, show_overflows = True
            )
            # Get the contribution from the FR method
            mc_bkg_fr = plotter.get_histogram(
                to_correct,
                '/mmt/final/Ntuple:' + selection + var + '_bkg_fr',
                rebin = rebin, show_overflows = True
            )
            mc_correct = mc_final - mc_bkg_fr
            corrected_mc_histos.append(mc_correct)

        signal = plotter.get_histogram(
            'VH120',
            '/mmt/final/Ntuple:' + selection + var + '_fin',
            rebin = rebin, show_overflows = True
        )

        stack = ROOT.THStack(selection + var + "FR_FINAL", "Final #mu#mu#tau selection")
        for histo in corrected_mc_histos:
            stack.Add(histo.th1, 'hist')
            stack.Add(histo.th1, 'hist')

        styling.apply_style(data_bkg_fr, **samplestyles.SAMPLE_STYLES['ztt'])
        stack.Add(data_bkg_fr.th1, 'hist')

        signal = signal*5
        stack.Add(signal.th1, 'hist')
        stack.Draw()
        data.Draw('pe,same')
        stack.SetMaximum(max(stack.GetMaximum(), data.GetMaximum())*1.5)
        stack.GetXaxis().SetTitle(x_title)

        saveplot(selection + var + '_fr')
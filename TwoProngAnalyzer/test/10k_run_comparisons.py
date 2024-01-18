import ROOT
import sys
import time

start_time = time.time()

xpix = 1600
ypix = 900
savedir = '10k_run/'
fileformat = '.png'

old_file = sys.argv[1] 
new_file = sys.argv[2]
vars_file = sys.argv[3]

old = []
new = []
dists = []

treename = "twoprongNtuplizer/fTree"
fo = open(old_file)
for line in fo:
  entry = {}
  chain = ROOT.TChain(treename)
  l = line.split()
  chain.Add(l[0])
  entry['tag'] = l[1]
  entry['tree'] = chain
  old.append(entry)

print "old"
print old

fn = open(new_file)
for line in fn:
  entry = {}
  chain = ROOT.TChain(treename)
  l = line.split()
  chain.Add(l[0])
  entry['tag'] = l[1]
  entry['tree'] = chain
  new.append(entry)

print "new"
print new

fv = open(vars_file)
for line in fv:
  dist = {}
  l = line.split()
  if len(l) == 0: continue
  if len(l) == 4:
    dist['var'] = l[0]
    dist['nbins'] = l[1]
    dist['low'] = l[2]
    dist['high'] = l[3]
    dist['cut'] = ""
  elif len(l) >= 5:
    dist['var'] = l[0]
    dist['nbins'] = l[1]
    dist['low'] = l[2]
    dist['high'] = l[3]
    cut_string = ""
    ll = l[4:len(l)]
    for p in ll:
      cut_string += p
    dist['cut'] = cut_string
  dists.append(dist)

print "dists"
for dist in dists:
  print dist

print "\nnow do plotting"
for i in range(len(old)):
  old_tree = old[i]['tree']
  new_tree = new[i]['tree']
  for dist in dists:
    print dist
    var = dist['var']
    nbins = int(dist['nbins'])
    low = float(dist['low'])
    high = float(dist['high'])
    histtitle = var + '_' + old[i]['tag'] + '_' + new[i]['tag']
    hist1 = ROOT.TH1D('hist1',histtitle,nbins,low,high)
    hist2 = ROOT.TH1D('hist2',histtitle,nbins,low,high)
    draw_string1 = var + " >> hist1" 
    draw_string2 = var + " >> hist2" 
    cut_string = dist['cut']
    print draw_string1, cut_string
    print draw_string2, cut_string
    old_tree.Draw(draw_string1, cut_string, "goff")
    new_tree.Draw(draw_string2, cut_string, "goff")
    c = ROOT.TCanvas('c1','c1',xpix,ypix)
    hist1.SetLineColor(ROOT.kBlack)
    hist2.SetLineColor(ROOT.kRed)
    hist1.SetStats(0)
    hist2.SetStats(0)
    hist1.Draw()
    hist2.Draw('same')
    leg = ROOT.TLegend(0.65,0.75,0.9,0.9)
    leg.AddEntry(hist1,'old','l')
    leg.AddEntry(hist2,'new','l')
    leg.Draw('same')
    savename = savedir + var + '_' + old[i]['tag'] + '_' + new[i]['tag'] + fileformat
    c.SaveAs(savename)

end_time = time.time()

print end_time - start_time

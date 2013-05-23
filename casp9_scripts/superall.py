from pymol import cmd

def sa(intra=False,rainbow=True):
  """
  Superimpose all open models onto the first one.  This may not work well with selections.  Option intra can be set to True to enable intra_fit first, for working with multi-state (nmr) pdbs.
  """
  AllObj=cmd.get_names("all")
  for x in AllObj:
    print(AllObj[0],x)
    if intra==True:
      cmd.intra_fit(x)
    if rainbow==True:
      cmd.util.chainbow(x)
    cmd.align(x,AllObj[0])
    cmd.zoom()

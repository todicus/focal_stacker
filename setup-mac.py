from setuptools import setup

setup(
  app=["focal_stacker_GUI.py"],
  options={"py2app":
    {"argv_emulation": True, "includes": ["sip", "PyQt4._qt"]}
  },
  setup_requires=["py2app"]) 
  
#from setuptools import setup
#
#setup(
#  app=["liftr/liftr.py"],
#  options={"py2app":
#    {"argv_emulation": True, "includes": ["sip", "PyQt4._qt"]}
#  },
#  setup_requires=["py2app"]) 
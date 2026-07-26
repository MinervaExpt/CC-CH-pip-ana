import os, tarfile, optparse, shutil, subprocess, errno, glob
import datetime as dt
import os.path




###############################################################################
# Constants/Default Args
###############################################################################
# Scripts, Files, and Dirs
kGRID_SCRIPT = os.getenv("PWD") + "/grid_ccpi_macro.sh"
kTOPDIR = os.getenv("TOPDIR")
kANATUPLE_DIR = "/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/"
#kANATUPLE_DIR = "/pnfs/minerva/persistent/DataPreservation/p8/FullDetector/"
#kANATUPLE_DIR = "/pnfs/minerva/persistent/users/zdar/"
kOUTDIR = "/pnfs/{EXPERIMENT}/scratch/users/{USER}/TestMAD/".format(
    EXPERIMENT=os.getenv("EXPERIMENT"), USER=os.getenv("USER")
)
kCACHE_PNFS_AREA = "/pnfs/{EXPERIMENT}/scratch/users/{USER}/grid_cache/".format(
    EXPERIMENT=os.getenv("EXPERIMENT"), USER=os.getenv("USER")
)
kTARBALL_LOCATION = "/pnfs/{EXPERIMENT}/scratch/users/{USER}/tarballs/".format(
    EXPERIMENT=os.getenv("EXPERIMENT"), USER=os.getenv("USER")
)
kEV_SEL_MACRO = "event_selection/runEventSelectionGrid.C+"
kMC_INPUTS_MACRO = "xsec/makeCrossSectionMCInputs.C+"
kBACKGROUND_BREAKDOWN = "studies/runBackgrounds.C+"
# Grid Stuff
kMINERVA_RELEASE = os.getenv("MINERVA_RELEASE")
kMEMORY = "10GB"
kGRID_OPTIONS = (
    "--group minerva "
    "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC "
#    "--lines='+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\\\"' "
    "--role=Analysis "
#    "--singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest "
#    "--OS=SL7 " # change to SL7 when submitting from sl7 machines.
)

# Misc
kPLAYLISTS = [
    "me1A",
    "me1B",
    "me1C",
    "me1D",
    "me1E",
    "me1F",
    "me1G",
    "me1L",
    "me1M",
    "me1N",
    "me1O",
    "me1P",
]
#kPLAYLISTS = [
#    "1A",
#    "1B",
#    "1C",
#    "1D",
#    "1E",
#    "1F",
#    "1G",
#    "1L",
#    "1M",
#    "1N",
#    "1O",
#    "1P",
#]
kFILETAG = ""


###############################################################################
# Helper Functions
###############################################################################
# mkdir -p
def MakeDirectory(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def CopyFile(source, destination):
    destination_full_path = destination + "/" + source
    if not os.path.exists(destination_full_path):
        shutil.copy(source, destination_full_path)
        print("copying " + source + " --> " + destination_full_path)
        return destination_full_path


def IFDHMove(source, destination):
    cmd = "mv " + source + " " + destination
    status = subprocess.call(cmd, shell=True)
    destination_full_path = destination + "/" + source
    return destination_full_path


def IFDHCopy(source, destination):
   # cmd = "ifdh cp " + source + " " + destination + "/" + source
    cmd = "cp " + source + " " + destination + "/" + source
    status = subprocess.call(cmd, shell=True)
    destination_full_path = destination + "/" + source
    return destination_full_path


# Tar up the given source directory.
# Right now, we only need Ana/ so skip everything else.
def MakeTarfile(source_dir, tag):
    tarfile_name = "granados_" + tag + ".tar.gz"
    tarfile_Name = "granados_" + tag + ".gz"
    # Do it
    tar = tarfile.open(tarfile_name, "w:gz")
    for i in os.listdir(source_dir):
        print(i)
        if (
            i == "Rec"
            or i == "Tools"
            or i == "Personal"
           # or i == "GENIEXSecExtract"
            or i == "Systemtics"
            or ("*.root" in i)
            or ("*.png" in i)
            or (".git" in i)
            or (".o" in i)
            or (".so" in i)
            or (".d" in i)
            or (".pcm" in i)
            or ("tar.gz" in i)
        ):
            continue
        print(source_dir + "/" + i)
        tar.add(source_dir + "/" + i, i)
    tar.close()

    # It is done. Send it to scratch.
    tarfile_fullpath = IFDHMove(tarfile_name, kTARBALL_LOCATION)

    return tarfile_name, tarfile_fullpath


def MakeUniqueProcessingID(tag):
    processing_id = "{TAG}{DAY}-{TIME}".format(
        TAG=tag, DAY=dt.date.today(), TIME=dt.datetime.today().strftime("%H%M%S")
    )
    return processing_id


def GetOptions():
    parser = optparse.OptionParser(usage="usage: %prog [options]")
    grid_group = optparse.OptionGroup(parser, "Grid Options")

    # grid args
    grid_group.add_option("--out_dir", default=kOUTDIR, help="Default = %default.")
    grid_group.add_option("--filetag", default=kFILETAG)
    grid_group.add_option("--tarfile", default="")
    grid_group.add_option("--memory", default=kMEMORY)
    grid_group.add_option("--ev_sel", action="store_true")
    grid_group.add_option("--mc_xsec_inputs", action="store_true")
    grid_group.add_option("--bg_breakdown", action="store_true")

    # job args
    job_group = optparse.OptionGroup(parser, "Job Options")

    job_group.add_option(
        "--no_truth",
        action="store_false",
        dest="do_truth",
        default=True,  # default is to DO truth
        help="Don't run over truth. Default: DO run over truth.",
    )
    job_group.add_option(
        "--no_systs",
        action="store_false",
        dest="do_full_systematics",
        default=True,  # default: DO systs
        help="Don't do full systematics. Default: DO systematics",
    )
    job_group.add_option(
        "--test_playlist",
        action="store_true",
        dest="do_test_playlist",
        default=False,  # default is Not use test playlist
        help="It use a test playlist. Default: It doesn't use the test playlist.",
    )
    job_group.add_option(
        "--signal_definition",
        default=0,
        help="0 = OnePiTracked,1 = OnePiTracked+Untracked, 2 = OnePiNoW Tracked+Untracked, 3 = NPi Tracked+Untracked, 4 = NPiNoW Tracked+Untracked, 5 = Nuke Tracked, 6 = OnePiTracked+Untracked for thetapi, 7 = = OnePiUntracked",
    )
    job_group.add_option(
        "--playlists",
        default="ALL",
        help="ALL, or ME1A, etc. Maybe a list someday but not now.",
    )
    job_group.add_option("--run", default=[], help="Specify a specific run number")

    parser.add_option_group(grid_group)
    parser.add_option_group(job_group)

    options, remainder = parser.parse_args()

    # require a macro
    if (
        options.ev_sel == options.mc_xsec_inputs
        and options.mc_xsec_inputs == options.bg_breakdown
    ):
        print("Pick a macro!")
        quit()
    elif options.ev_sel:
        options.macro = kEV_SEL_MACRO
    elif options.mc_xsec_inputs:
        options.macro = kMC_INPUTS_MACRO
    elif options.bg_breakdown:
        options.macro = kBACKGROUND_BREAKDOWN
    else:
        pass

    # fix file tag underscore
    if options.filetag != kFILETAG and not options.filetag[:1] == "_":
        options.filetag = "_" + options.filetag

    return options


###############################################################################
# Main
###############################################################################
def main():
    options = GetOptions()

    # A unique string for this job
    processing_id = MakeUniqueProcessingID(options.filetag)

    # Make out_dir
    print("Outdir (top) is " + options.out_dir)
    out_dir = options.out_dir + "/" + processing_id
    MakeDirectory(out_dir)

    # Make tarfile and pass to resilient
    if options.tarfile:
        tarfile = options.tarfile.split("/")[-1]
        tarfile_fullpath = options.tarfile
    else:
        tarfile, tarfile_fullpath = MakeTarfile(kTOPDIR, processing_id)

    print("\nUsing tarfile: " + tarfile_fullpath)

    # Let's send the grid script to pnfs first:
    cache = kCACHE_PNFS_AREA + "/" + processing_id
    print("sending grid macro to " + cache)
    MakeDirectory(cache)
    grid_script = IFDHCopy("grid_ccpi_macro.sh", cache)

    if options.run:
        print("\nSubmitting run: " + options.run)

    # Loop playlists, anatuples, and submit
    for i_playlist in kPLAYLISTS:
        do_this_playlist = i_playlist == options.playlists or options.playlists == "ALL"
        if not do_this_playlist:
            continue

        print("Using tuples from" + kANATUPLE_DIR)

        # loop anatuples
#        list_of_anatuples = glob.glob(
#            kANATUPLE_DIR + "/Playlist{0}/*".format(i_playlist)
#        )
        list_of_anatuples = glob.glob(
            kANATUPLE_DIR + "/Merged_mc_ana_{0}_DualVertex_p4/*".format(i_playlist)
        )
#        list_of_anatuples = glob.glob(
#            kANATUPLE_DIR + "/Merged_mc_ana_{0}_DualVertex_preP6/*".format(i_playlist)
#        )
        #        list_of_anatuples = glob.glob(kANATUPLE_DIR + "/Merged_mc_ana_{0}_DualVertex_preP4_13June2023/*".format(i_playlist))
        for anatuple in list_of_anatuples:
            if not ("MasterAnaDev" in anatuple) or not (".root" in anatuple):
                continue

            run = anatuple[-22:-14]
            # run = anatuple[-13:-5] # merging with outdated custom method
            run = run.lstrip("0")
            if options.run and (run not in options.run):
                continue
            print(anatuple)
            print(run)

            def XROOTDify(anatuple):
#                return anatuple.replace(
#                    "/pnfs/", "root://fndcadoor.fnal.gov/pnfs/fnal.gov/usr/"
#                )
                return anatuple.replace(
                    "/pnfs/", "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/"
                )

            anatuple = XROOTDify(anatuple)
            if options.mc_xsec_inputs:
                macro = options.macro
                macro += (
                    '({SIGNAL_DEFINITION},\\\\\\"{PLAYLIST}\\\\\\",{DO_FULL_SYST},'
                    '{DO_TRUTH},{DO_TEST},{DO_GRID},\\\\\\"{TUPLE}\\\\\\",{RUN})'.format(
                        SIGNAL_DEFINITION=options.signal_definition,
                        PLAYLIST=i_playlist,
                        DO_FULL_SYST="true" if options.do_full_systematics else "false",
                        DO_TRUTH="true" if options.do_truth else "false",
                        DO_TEST="false" if options.do_test_playlist else "true",
                        DO_GRID="true",
                        TUPLE=anatuple,
                        RUN=run,
                    )
                )

            if options.bg_breakdown:
                macro = options.macro
                macro += (
                    '({SIGNAL_DEFINITION},\\\\\\"{PLAYLIST}\\\\\\",'
                    '{DO_GRID},\\\\\\"{TUPLE}\\\\\\",{RUN})'.format(
                        SIGNAL_DEFINITION=options.signal_definition,
                        PLAYLIST=i_playlist,
                        DO_GRID="true",
                        TUPLE=anatuple,
                        RUN=run,
                    )
                )

            macro = '"' + macro + '"'

            print(macro)

            # Prepare Submit Command
            submit_command = (
                "jobsub_submit {GRID} --memory {MEMORY} --expected-lifetime=24h "
                "-d OUT {OUTDIR} "
                "-L {LOGFILE} "
                "-e MACRO={MACRO} "
                "-e TARFILE={TARFILE} "
                "-f dropbox://{TARFILE_FULLPATH} "
                #                "--tar_file_name dropbox://{TARFILE_FULLPATH} "
                "--use-pnfs-dropbox "
                "file://{GRID_SCRIPT}".format(
                    GRID=kGRID_OPTIONS,
                    MEMORY=options.memory,
                    OUTDIR=out_dir,
                    LOGFILE=out_dir + "/log{0}.txt".format(run),
                    MACRO=macro,
                    TARFILE=tarfile,
                    TARFILE_FULLPATH=tarfile_fullpath,
                    GRID_SCRIPT=grid_script,
                )
            )  # submit_command

            # Ship it
            print("\nSubmitting to grid:\n" + submit_command + "\n")
            status = subprocess.call(submit_command, shell=True)


if __name__ == "__main__":
    main()

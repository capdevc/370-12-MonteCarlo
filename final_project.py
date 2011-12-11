#!/usr/bin/python

# Assignment 12: Monte Carlo
# Cristian Capdevila 11/27/11
# See README

from numpy import *
from numpy.random import *
from pylab import *
import multiprocessing
import os, sys
import ConfigParser
import datetime


# Globals for callbacks
infected = array([0], dtype=float64)
vaccinated = array([0], dtype=float64)
recovered = array([0], dtype=float64)
infhist = []
datfile = 0


def sim(tau, nu, size, patients, k, save, runnum):
    # The main simulation function
    # save is set to 1 if this run is being saved for animation.
    room = zeros((size, size))
    room[randint(size)][randint(size)] = 1  # place infected
    inf = [1]  # local lists of daily totals
    vac = [0]
    rec = [0]
    infctr = 1  # local counters of each status
    vacctr = 0
    recctr = 0
    maps = []  # list of daily room maps

    # this randomly places any empty beds
    for x in range(size**2 - patients):
        slot_found = 0
        while slot_found == 0:
            i = randint(size)
            j = randint(size)
            if room[i][j] == 0:
                room[i][j] = -3
                slot_found = 1
    newroom = room.copy()
    if save: maps.append(room)

    # main loop that runs through patients as long as at least one is infected
    # -3 = empty bed
    # -2 = vaccinated
    # -1 = recovered
    #  0 = susceptible
    # >0 = days infected
    while infctr > 0:
        for i in range(0, size):
            for j in range(0, size):
                if (room[i][j] >= k): 
                    newroom[i][j] = -1
                    infctr -= 1
                    recctr += 1
                elif room[i][j] > 0: newroom[i][j] += 1
                elif room[i][j] == 0: 
                    if rand() < nu: 
                        newroom[i][j] = -2
                        vacctr += 1
                    elif infect(room, i, j, tau, size): 
                        newroom[i][j] = 1
                        infctr += 1
                else: newroom[i][j] = room[i][j]
        room = newroom.copy() # update the map all at once
        if save: maps.append(room)
        inf.append(infctr)
        vac.append(vacctr)
        rec.append(recctr)
    if save: return ((inf, vac, rec, runnum), maps)
    else: return (inf, vac, rec, runnum) # return lists of daily counts

def infect(room, i, j, tau, size):
    t = 1  # t is 0 if patient infected, 1 otherwise
    if i > 0:
        if room[i-1, j] > 0: t *= (rand() > tau) # north
    if i < size-1:
        if room[i+1, j] > 0: t *= (rand() > tau) # south
    if j > 0:
        if room[i, j-1] > 0: t *= (rand() > tau) # west
    if j < size-1:
        if room[i, j+1] > 0: t *= (rand() > tau) # east
    return 1 if t == 0 else 0 # return 1 if infected

def cb(counts):
    # callback function for processes
    # results are aggregated into the global lists here
    global infected, vaccinated, recovered, infhist
    delta = len(counts[0])- len(infected)
    # counts is the tuple returned from sim
    inf = array(counts[0])
    vac = array(counts[1])
    rec = array(counts[2])
    runnum = counts[3]
    if datfile:  # print to file
        datfile.write('Run ' + str(runnum) + ' inf counts: ' + str(inf) + '\n')
        datfile.write('Run ' + str(runnum) + ' vac counts: ' + str(vac) + '\n')
        datfile.write('Run ' + str(runnum) + ' rec counts: ' + str(rec) + '\n')
    if delta > 0: # global list is shorter than this run
        infected.resize((len(inf)))
        # carry over last day's vaccinated and recovered counts io future days
        vaccinated = concatenate((vaccinated, vaccinated[-1] * ones(delta)))
        recovered = concatenate((recovered, recovered[-1] * ones(delta)))
    elif delta < 0: # this run is shorter than global list
        inf.resize((len(infected)))
        vac = concatenate((vac, vac[-1] * ones(-delta)))
        rec = concatenate((rec, rec[-1] * ones(-delta)))
    # increment the global lists with info from last run and add the infected
    # details to a history list to do statistics on. The callbacks all run in
    # the main process, so no need to worry about data races on globals
    infected += inf
    vaccinated += vac
    recovered += rec
    infhist.append(inf.tolist())

def colorize(maps):
    # takes a list of room maps and converts the infected states into rgb
    images = []
    for map in maps:
        image = zeros((size, size, 3), dtype=float64)
        i = 0
        for row in map:
            j = 0
            for x in row:  # substitute RGB values
                if (x == -3): image[i][j] = array([0.0, 0.0, 0.0])
                elif (x == -2): image[i][j] = array([0.0, 0.0, 1.0])
                elif (x == -1): image[i][j] = array([0.0, 1.0, 0.0])
                elif (x == 0): image[i][j] = array([0.75, 0.75, 0.75])
                else: image[i][j] = array([1.0, (1.0*x)/k, 0.0])
                j += 1
            i += 1
        images.append(image)
    return images

def imagegen(maps):
    # generate the images and movie from the saved run
    images = colorize(maps)
    # call imshow and save png files
    files =[]
    fig = figure(figsize=(7,7))
    ax = fig.add_subplot(111)
    i = 1
    for image in images:
        ax.cla()
        ax.imshow(image, cmap='spectral', interpolation='nearest')
        fname = 'mcarlo%03d.png' %i
        fig.savefig(fname)
        files.append(fname)
        i += 1
    # use mencoder to make a movie, about 10 seconds long
    fps = len(maps)/10
    print "Generating mcarlo.avi..."
    os.system("mencoder 'mf://mcarlo*.png' -mf type=png:fps=%i -ovc lavc " \
        "-lavcopts vcodec=wmv2 -oac copy -o mcarlo.avi > /dev/null 2>&1" % fps)

# main process
if __name__ == '__main__':
    # Read Parameters from config file
    config = ConfigParser.RawConfigParser()
    config.read('final_project.conf')
    tau = config.getfloat('Parameters', 'Tau')
    nu = config.getfloat('Parameters', 'Nu')
    size = config.getint('Parameters', 'Room Side Length')
    patients = config.getint('Parameters', 'Patient Count')
    k = config.getint('Parameters', 'Infection Duration')
    runs = config.getint('Parameters', 'Sim Runs')
    savedrun = config.getint('Parameters', 'Animated Run')
    processes = config.getint('Parameters', 'Processes')
    printdata = config.getint('Parameters', 'Print')

    # bounds check the parameters
    if tau < 0 or tau > 1:
       print "Tau:%0.3f. Valid settings are between 0.0 and 1.0." % tau
       sys.exit()
    if nu < 0 or nu > 1:
       print "Nu:%0.3f. Valid settings are between 0.0 and 1.0." % nu
       sys.exit()
    if size < 1:
        print "Room Side Length:%d. Must be greater than 0." % size
        sys.exit()
    if patients < 0:
        print "Patient Count:%d. Must be greater than 1." % patients
        sys.exit()
    if patients > size **2: 
        print "Patient Count: %d. Setting to max for room size, %d." \
                % (patients, size**2)
        patients = size**2
    if patients == 0: patients = size **2
    if k < 0:
        print "Infection Duration:%d. Must be greater than 0." % k
        sys.exit()
    if runs < 1:
        print "Sim Runs:%d. Must be greater than 0." % runs
        sys.exit()
    if savedrun > runs:
        print "Animated Run:%d. Must be less that Sim Runs. Selecting a run " \
                "to animate at random." % savedrun
        savedrun = -1
    if savedrun == -1 or savedrun >= runs:
        savedrun = randint(1, runs+1)
        print "Selected run", savedrun, "to animate."
    if processes < 0:
        print "Processes:%d. Must be greater than 0. Setting to number of "\
               "cores available." % processes
        processes = multiprocessing.cpu_count()
    savedmaps = multiprocessing.Queue() # queue holding the saved run

    if processes == 0: processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes)  # set up a process pool

    datfile = 0
    if printdata:
        datfile = open('mcarlo.txt', 'w')
        print "Saving raw data to mcarlo.txt."
    print "Running %i simulations with %d patients in a %dx%d room in %d " \
            "processes..." % (runs, patients, size, size, processes)
    time1 = datetime.datetime.now() # start a timer
    result = []
    for runnum in range(1, runs+1):
        # here we assign a sim job to the process pool
        if runnum == savedrun: result = pool.apply_async(sim, (tau, nu, size, 
            patients, k, 1, runnum)) # save this run
        else: pool.apply_async(sim, (tau, nu, size, patients, k, 0, runnum),
            callback=cb)
    pool.close()
    pool.join()  # join the processes and finish the callbacks
    print "Elapsed time in processes: ", datetime.datetime.now() - time1

    if savedrun:
        # Add the saved run to stats by calling cb and then generate images
        print "Generating images..."
        tup, maps = result.get()
        cb(tup)
        imagegen(maps)

    # get stats on infection duration
    days = []
    for run in infhist:
        days.append(len(run))
    maxdays = max(days)

    # square up the list of infected data and calculate std deviation
    # to plot as error bars
    for list in infhist:
        list.extend([0] * (maxdays - len(list)))
    infhista = array(infhist)
    errors = []
    for i in range(maxdays):
        errors.append(std(infhista[0:,i]))

    # Convert the data to averages
    infected = infected/runs
    vaccinated = vaccinated/runs
    recovered = recovered/runs

    if printdata:  # print the data we're plotting
        print "\nAverage daily infected:", infected
        print "Average daily vaccinated", vaccinated
        print "Average daily recovered", recovered

    # Print the most interesting data
    print "\nShortest run: ", min(days), "days."
    print "Longest run:  ", maxdays, "days."
    print "Median run:   ", int(median(days)), "days."
    print "Mean of runs: ", mean(days), "days."
    print "Std deviation: %0.3f days." % std(days)

    # plot the graphs last since show() blocks
    xdat = range(len(infected))
    clf()
    plot(xdat, vaccinated, 'b-', label='Avg. # Vaccinated', lw=1.0)
    plot(xdat, recovered, 'g-', label='Avg. # Recovered', lw=1.0)
    errorbar(xdat, infected, yerr=errors, ecolor='r', color='r', 
            label='Avg. # Infected w/ StdDev', ls='-', lw=1.0)
    legend(loc=2)
    title("%i runs with Tau: %0.3f and Nu: %0.3f" % (runs, tau, nu),
            fontsize='x-large')
    xlabel("Elapsed Time (Days)")
    ylabel("Patients")
    print "Finished."
    show()

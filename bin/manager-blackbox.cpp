/*
 * manager-external.cpp
 *
 * Created on: 21/06/2016
 *
 * This executable starts a manager which can monitor applications not written
 * with the Nornir framework.
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of nornir.
 *
 *  nornir is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  nornir is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with nornir.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <sys/mman.h>
#include <sys/wait.h>
#include <stdlib.h>
#include "../src/manager.hpp"
#include "../src/manager-multi.hpp"
#include <tclap/CmdLine.h>

using namespace nornir;

#undef DEBUG
#undef DEBUGB

#define DEBUG_MBB

#ifdef DEBUG_MBB
#define DEBUG(x) do { cerr << "[BlackBox Manager] " << x << endl; } while (0)
#define DEBUGB(x) do {x;} while(0)
#else
#define DEBUG(x)
#define DEBUGB(x)
#endif

static volatile bool *started, *handlerCreated;

struct ScheduledProgram;

typedef struct ScheduledProgram{
    friend std::ostream& operator <<(std::ostream& os, const ScheduledProgram& instance);
    double start;
    std::string parametersFile;
    std::vector<std::string> program;
    pid_t pid;
}ScheduledProgram;

std::ostream& operator <<(std::ostream& os, const ScheduledProgram& instance){
    os << "Starts at: [" << instance.start << "] ";
    os << "Parameters: [" << instance.parametersFile << "] ";
    os << "Program: [";
    for(size_t i = 0; i < instance.program.size(); i++){
        os << instance.program[i] << " ";
    }
    os << "]";
    return os;
}

static std::string pidToString(pid_t pid){
    stringstream out;
    out << pid;
    return out.str();
}

int main(int argc, char * argv[]){
    std::vector<long> pids;
    std::vector<std::string> pidsParameters;
    std::vector<ScheduledProgram> scheduledPrograms;
    try {
        TCLAP::CmdLine cmd("Runs Nornir on already existing processes, coordinating "
                           "their requirements. It is also possible to specify new "
                           "programs to be ran. One parameter between --pid or "
                           "--schedule must be specified.", ' ', "1.0");
        TCLAP::MultiArg<long> pidFlag("p", "pid", "PID of a running process.", false, "long", cmd);
        TCLAP::MultiArg<std::string> parametersFlag("r", "parameters", "Parameters for a running process.", false, "string", cmd);
        TCLAP::ValueArg<std::string> scheduleFlag("s", "schedule", "Name of the file containing the schedule of "
                                                  "the programs to be run. The file has multiple lines, each one with "
                                                  "the following syntax:\n"
                                                  "StartTime ParametersFile ProgramPath ProgramArguments\n"
                                                  "The start time is a relative offset starting from the start of this "
                                                  "executable. The lines must be sorted increasly by start time.", 
                                                  false, "", "string", cmd);
        cmd.parse(argc, argv);

#if 1
        if(pidFlag.getValue().size()){
            std::cerr << "[ERROR] -p not still supported (due to problems in reading performance counters for non-child processes." << std::endl;
            return -1;
        }
#endif

        /* Validate parameters. */
        if(!pidFlag.getValue().size() &&
           !scheduleFlag.getValue().compare("")){
            std::cerr << "[ERROR] One between --run, --pid or --schedule must be specified." << std::endl;
            return -1;
        }

        if(pidFlag.getValue().size() !=
           parametersFlag.getValue().size()){
            std::cerr << "[ERROR] Number of pids and number of specified parameters must be equal." << std::endl;
            return -1;
        }

        /* Load schedule. */
        pids = pidFlag.getValue();
        pidsParameters = parametersFlag.getValue();
        std::ifstream scheduleFile(scheduleFlag.getValue());
        std::string line;
        while(std::getline(scheduleFile, line)){
            std::vector<std::string> fields;
            fields = mammut::utils::split(line, ' ');
            if(fields.size() < 3){
                std::cerr << "Invalid scheduled line: " << line << ". You must at least "
                             "specify the scheduled time, parameters file and program name." << std::endl;
                return -1;
            }
            ScheduledProgram sp;
            sp.start = mammut::utils::stringToDouble(fields[0]);
            sp.parametersFile = fields[1];
            sp.program = std::vector<std::string>(fields.begin() + 2, fields.end());
            DEBUG("Read schedule: " << sp);
            scheduledPrograms.push_back(sp);
        }
        scheduleFile.close();
    }catch (TCLAP::ArgException &e){
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        return -1;
    }

    ManagerMulti mm;
    bool multiManagerNeeded = false;
    if(pids.size() + scheduledPrograms.size() > 1){
        multiManagerNeeded = true;
        mm.start();
        DEBUG("Multi manager started.");
    }

    /* First we add the already running pid. */
    for(size_t i = 0; i < pids.size(); i++){
        Parameters p(pidsParameters.at(i));
        std::string logPrefix = pidToString(pids.at(i));
        Observer o(logPrefix + "_stats.csv",
                   logPrefix + "_calibration.csv",
                   logPrefix + "_summary.csv");
        p.observer = &o;
        ManagerBlackBox* m = new ManagerBlackBox(pids.at(i), p);
        if(multiManagerNeeded){
            mm.addManager(m);
            DEBUG("Added PID: " << pids.at(i) << " to the multimanager.");
        }else{
            m->start();
            DEBUG("Started PID: " << pids.at(i));
            m->join();
            delete m;
        }
    }

    /* Then we start following the schedule. */
    started = (volatile bool*) mmap(NULL, sizeof(volatile bool), PROT_READ | PROT_WRITE,
                           MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    handlerCreated = (volatile bool*) mmap(NULL, sizeof(volatile bool), PROT_READ | PROT_WRITE,
                                  MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    *started = false;
    *handlerCreated = false;
    double lastStart = 0;

    for(size_t i = 0; i < scheduledPrograms.size(); i++){
        ScheduledProgram sp = scheduledPrograms[i];
        sleep(sp.start - lastStart);
        lastStart = sp.start;

        void* mmem = (void*) mmap(NULL, sizeof(ManagerBlackBox), PROT_READ | PROT_WRITE,
                                  MAP_SHARED | MAP_ANONYMOUS, -1, 0);

        pid_t pid = fork();
        if(pid == -1){
            std::cerr << "Fork failed." << std::endl;
            return -1;
        }else if(pid){
            /* Father - Manager. */
            sp.pid = pid;
            while(!*started){;}
            Parameters p(sp.parametersFile);
            stringstream out;
            out << sp.start;
            std::string logPrefix = out.str() + "_" + mammut::utils::split(sp.program.at(0), '/').back();
            Observer o(logPrefix + "_stats.csv",
                       logPrefix + "_calibration.csv",
                       logPrefix + "_summary.csv");
            p.observer = &o;
            while(!*started){;}
            ManagerBlackBox* m = new (mmem) ManagerBlackBox(p.mammut.getInstanceTask()->getProcessHandler(pid), p);
            *handlerCreated = true;
            if(multiManagerNeeded){
                mm.addManager(m);
                // Wait for the child to reset the flags.
                while(*handlerCreated){;}
            }else{
                m->start();
                m->join();
                waitpid(pid, NULL, 0);
                m->~ManagerBlackBox(); // Because created with placement new
            }
            DEBUG("Started scheduled: " << sp.program[0] << " with pid " << pid);
        }else{
            /* Child - Application */
            *started = true;
            // If the handler is not created, we would not catch the counters of
            // the threads/processes created by this process
            while(!*handlerCreated){;}
            extern char** environ;
            char** arguments = new char*[sp.program.size() + 1];
            for(size_t i = 0; i < sp.program.size(); i++){
                arguments[i] = new char[sp.program[i].size() + 1];
                strcpy(arguments[i], sp.program[i].c_str());
            }
            arguments[sp.program.size()] = NULL;
            DEBUG("Running program " << sp.program[0]);
            /* Resets flag for next scheduled application. */
            *started = false;
            *handlerCreated = false;
            if(execve(arguments[0], arguments, environ) == -1){
                std::cerr << "Impossible to run the specified executable. Terminating the manager." << std::endl;
                ((ManagerBlackBox*) mmem)->terminate();
                return -1;
            }
            return -1; // execve never returns (except in the error case above).
        }
    }

    DEBUG("All specified programs have been added to the monitoring system. Waiting for their termination...");
    /* All programs scheduled, now wait for termination and clean. */
    if(multiManagerNeeded){
        for(size_t i = 0; i < pids.size() + scheduledPrograms.size(); i++){
            Manager* m;
            while((m = mm.getTerminatedManager()) == NULL){
                sleep(1);
            }
            DEBUG("One program terminated.");
            // TODO:
            //if(scheduledProgram(pid)){
            //    m->~ManagerBlackBox(); // Because we used placement new
            //}else{
            //    delete m;
            //}
        }
    
        for(size_t i = 0; i < scheduledPrograms.size(); i++){
            waitpid(scheduledPrograms[i].pid, NULL, 0);
        }
    }
}


/*
 * monitor.hpp
 *
 * Created on: 03/07/2016
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

#ifndef NORNIR_CLIENT_HPP_
#define NORNIR_CLIENT_HPP_

#include "external/orlog/src/orlog.hpp"
#include "external/orlog/src/external/cppnanomsg/nn.hpp"
#include "external/orlog/src/external/nanomsg/src/pair.h"

#define EXTERNAL_CHANNEL_NAME "ipc:///tmp/nornir.ipc"

namespace nornir{
class ExternalApplication{
private:
    nn::socket _channel;
    int _chid;
    orlog::Application* _app;
public:
    /**
     * Creates a client for interaction with a local server.
     * @param parameters The file containing the Nornir parameters.
     */
    ExternalApplication(const std::string& parametersFile);

    /**
     * Creates a client for interaction with a remote server.
     * @param parameters The file containing the Nornir parameters.
     * @param serverAddress The address of the remote Nornir server.
     * @param port The port where the remote Nornir server is listening.
     */
    ExternalApplication(const std::string& parametersFile, const std::string& serverAddress, uint port = 9999);

    /**
     * Destroys the handle.
     */
    ~ExternalApplication();

    /**
     * This function must be called at each loop iteration when the computation
     * part of the loop begins.
     */
    void begin();

    /**
     * This function must be called at each loop iteration when the computation
     * part of the loop ends.
     */
    void end();

    /**
     * This function must only be called once, when the parallel part
     * of the application terminates.
     */
    void terminate();
};

}

#endif /* NORNIR_CLIENT_HPP_ */

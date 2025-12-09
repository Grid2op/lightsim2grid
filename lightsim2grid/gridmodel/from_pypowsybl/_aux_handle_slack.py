# Copyright (c) 2025, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.


import numpy as np


def handle_slack_one_el(df_gen, gen_slack_id):
    if isinstance(gen_slack_id, (str, np.str_)):
        gen_slack_id_int = int((df_gen.index == gen_slack_id).nonzero()[0][0])
        gen_slack_weight = -float(df_gen.iloc[gen_slack_id_int]["p"])
    elif isinstance(gen_slack_id, tuple):
        if len(gen_slack_id) != 2:
            raise RuntimeError("When initializing gen_slack_id with a tuple, it should be "
                               "(gen_id_or_name, slack_weights). You provided a tuple of size "
                               f"{len(gen_slack_id)} != 2.")
        gen_slack_id_maybe_int, gen_slack_weight = gen_slack_id
        if not isinstance(gen_slack_id_maybe_int, int):
            gen_slack_id_int = handle_slack_one_el(df_gen, gen_slack_id_maybe_int)[0]
        else:
            gen_slack_id_int = gen_slack_id_maybe_int
        if float(gen_slack_weight) != gen_slack_weight:
            raise RuntimeError("When initializing the slack with a tuple, the second element of "
                               "each is the slack weight and should be convertible to a float.")
        gen_slack_weight = float(gen_slack_weight)
    else:
        try:
            gen_slack_id_int = int(gen_slack_id)
        except Exception as exc_:
            raise RuntimeError("'gen_slack_id' should be either an int or "
                                "a generator names") from exc_
        if gen_slack_id_int != gen_slack_id:
            raise RuntimeError("'gen_slack_id' should be either an int or a "
                                "generator names")
        gen_slack_id_int = gen_slack_id_int
        gen_slack_weight = 1.
    if not df_gen.iloc[gen_slack_id_int]["connected"] or abs(gen_slack_weight) < 1e-5:
        return None, None
    return gen_slack_id_int, gen_slack_weight


def handle_slack_iterable(df_gen, gen_slack_id):
    res_ids = []
    res_ws = []
    if isinstance(gen_slack_id, (list, tuple, set)):
        for el in gen_slack_id:
            tmp_id, tmp_w = handle_slack_one_el(df_gen, el)
            if tmp_id is None or tmp_w is None:
                # gen is disconnected
                continue
            res_ids.append(tmp_id)
            res_ws.append(tmp_w)
    elif isinstance(gen_slack_id, (dict)):
        for k, v in gen_slack_id.items():
            tmp_id, tmp_w = handle_slack_one_el(df_gen, (k, v))
            if tmp_id is None or tmp_w is None:
                # gen is disconnected
                continue
            res_ids.append(tmp_id)
            res_ws.append(tmp_w)
    else:
        raise RuntimeError("when intializing the slack with an iterable, make sure to "
                           "provide either a python list or a python dict.")
    return res_ids, res_ws
            
import pypowsybl as pypow
import pypowsybl.loadflow as pypow_lf


def get_same_slack(case_name):
    if case_name == "ieee9":
        return "VL1_0", 0
    if case_name == "ieee14":
        return "VL1_0", 0
    if case_name == "ieee30":
        return "VL1_0", 0
    if case_name == "ieee57":
        return "VL1_0", 0
    if case_name == "ieee118":
        return "VL69_0", 68
    if case_name == "ieee300":
        return "VL7049_0", 257
    
    raise RuntimeError(f"Unknown env {case_name}")


def get_pypowsybl_parameters(slack_voltage_level=None):    
    params = pypow_lf.Parameters(
        voltage_init_mode=pypow._pypowsybl.VoltageInitMode.UNIFORM_VALUES,
        transformer_voltage_control_on=False,
        use_reactive_limits=False,
        phase_shifter_regulation_on=False,
        twt_split_shunt_admittance=True,
        shunt_compensator_voltage_control_on=False,
        read_slack_bus=False,
        write_slack_bus=True,
        distributed_slack=False,
        dc_use_transformer_ratio=True,
        hvdc_ac_emulation=False,
        dc_power_factor=1.,
        provider_parameters={
            "useActiveLimits": "false",
            "useReactiveLimits": "false",
            "svcVoltageMonitoring": "false",
            "voltageRemoteControl": "false",
            "writeReferenceTerminals": "false",
            "voltagePerReactivePowerControl": "false",
            "generatorReactivePowerRemoteControl": "false",
            "secondaryVoltageControl": "false",
            # see https://github.com/powsybl/pypowsybl/issues/1127#issuecomment-3581713875
            'slackDistributionFailureBehavior' : 'LEAVE_ON_SLACK_BUS',
            # 'transformerVoltageControlMode' : 'WITH_GENERATOR_VOLTAGE_CONTROL',
            'plausibleActivePowerLimit': '5000'
            }
        )
    if slack_voltage_level is not None:
        params.provider_parameters["slackBusSelectionMode"] = "NAME"
        params.provider_parameters["slackBusesIds"] = f"{slack_voltage_level}"
    return params


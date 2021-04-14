function fmrwhy_util_createErrorMsg(function_name, identifier_txt, msg)

    errorStruct.identifier = [function_name ':' identifier_txt];
    errorStruct.message = msg;
    error(errorStruct);

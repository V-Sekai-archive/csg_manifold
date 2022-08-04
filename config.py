def can_build(env, platform):
    return not env["disable_3d"] and not env["platform"] == "android" and not env["platform"] == "javascript"


def configure(env):
    pass


def get_doc_classes():
    return [
        "CSGManifold3D",
    ]


def get_doc_path():
    return "doc_classes"

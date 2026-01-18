def classFactory(iface):
    from .terralink import TerraLinkPlugin
    return TerraLinkPlugin(iface)

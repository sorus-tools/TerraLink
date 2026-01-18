import os

from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction
from qgis.core import QgsApplication

from .processing_provider import TerraLinkProcessingProvider
try:
    from qgis.PyQt import sip
except Exception:  # pragma: no cover
    sip = None


class TerraLinkPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.plugin_dir = os.path.dirname(__file__)
        self.action = None
        self.provider = None

    def initGui(self):
        self.action = QAction("Run TerraLink", self.iface.mainWindow())
        icon_path = os.path.join(self.plugin_dir, "icon.png")
        icon = QIcon(icon_path)
        if not icon.isNull():
            self.action.setIcon(icon)
        self.action.triggered.connect(self.run)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("&TerraLink", self.action)

        if self.provider is None:
            self.provider = TerraLinkProcessingProvider(self)
            QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        if self.action is not None:
            self.iface.removeToolBarIcon(self.action)
            self.iface.removePluginMenu("&TerraLink", self.action)
        if self.provider is not None:
            try:
                registry = QgsApplication.processingRegistry()
                if registry is not None:
                    if sip is None or not sip.isdeleted(self.provider):
                        registry.removeProvider(self.provider)
            except RuntimeError:
                pass
            self.provider = None

    def run(self):
        from .terralink_dialog import TerraLinkDialog

        dialog = TerraLinkDialog(self.iface)
        dialog.exec_()

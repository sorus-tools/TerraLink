import os

from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction
from qgis.core import QgsApplication

from .processing_provider import TerraLinkProcessingProvider


class TerraLinkPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.plugin_dir = os.path.dirname(__file__)
        self.action = None
        self.provider = None

    def initGui(self):
        self.action = QAction("Run TerraLink (v1.1)", self.iface.mainWindow())
        icon_path = os.path.join(self.plugin_dir, "icon.png")
        icon = QIcon(icon_path)
        if not icon.isNull():
            self.action.setIcon(icon)
        self.action.triggered.connect(self.run)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("&TerraLink (v1.1)", self.action)

        if self.provider is None:
            self.provider = TerraLinkProcessingProvider(self)
            QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        if self.action is not None:
            self.iface.removeToolBarIcon(self.action)
            self.iface.removePluginMenu("&TerraLink (v1.1)", self.action)
        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None

    def run(self):
        from .linkscape_dialog import LinkscapeDialog

        dialog = LinkscapeDialog(self.iface)
        dialog.exec_()

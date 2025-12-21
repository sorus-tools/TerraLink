from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import QgsProcessingAlgorithm, QgsProcessingProvider


class TerraLinkLaunchAlgorithm(QgsProcessingAlgorithm):
    """Processing entry that opens the TerraLink dialog."""

    ALG_ID = "terralink_v1_1_corridor_generation"

    def __init__(self, plugin):
        super().__init__()
        self._plugin = plugin

    def name(self) -> str:
        return self.ALG_ID

    def displayName(self) -> str:
        return self.tr("TerraLink Corridor Generation")

    def group(self) -> str:
        # No subgroup: show directly under the SORUS provider.
        return ""

    def groupId(self) -> str:
        return ""

    def shortHelpString(self) -> str:
        return self.tr("Opens the TerraLink dialog so you can configure and run analyses.")

    def flags(self):
        return super().flags() | QgsProcessingAlgorithm.FlagNoThreading

    def tr(self, string: str) -> str:
        return QCoreApplication.translate("TerraLinkLaunchAlgorithm", string)

    def initAlgorithm(self, config=None):  # noqa: D401 - required override
        """No parameters to declare."""
        return

    def processAlgorithm(self, parameters, context, feedback):
        feedback.pushInfo(self.tr("Opening TerraLink dialog..."))
        self._plugin.run()
        return {}

    def createInstance(self):
        return TerraLinkLaunchAlgorithm(self._plugin)


class TerraLinkProcessingProvider(QgsProcessingProvider):
    """Processing provider that exposes the TerraLink launcher."""

    def __init__(self, plugin):
        super().__init__()
        self._plugin = plugin

    def id(self) -> str:
        return "terralink_v1_1"

    def name(self) -> str:
        return self.tr("TerraLink (v1.1)")

    def tr(self, string: str) -> str:
        return QCoreApplication.translate("TerraLinkProcessingProvider", string)

    def loadAlgorithms(self):
        self.addAlgorithm(TerraLinkLaunchAlgorithm(self._plugin))

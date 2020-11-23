// Copyright (c) Microsoft Corporation.
// Licensed under the MIT License.

#nullable enable

using System.Text;
using Microsoft.Jupyter.Core;
using Microsoft.Extensions.Logging;
using Newtonsoft.Json;
using static NumSharp.Slice;
using Microsoft.Quantum.IQSharp.Jupyter;
using System.Linq;
using System;

namespace Microsoft.Quantum.Diagnostics.Emulation
{
    public class DisplayableUnitaryOperatorToTextEncoder : IResultEncoder
    {
        private ILogger<DisplayableUnitaryOperatorToTextEncoder> logger;
        public string MimeType => MimeTypes.PlainText;

        public DisplayableUnitaryOperatorToTextEncoder(ILogger<DisplayableUnitaryOperatorToTextEncoder> logger)
        {
            this.logger = logger;
        }

        public EncodedData? Encode(object displayable)
        {
            if (displayable is DisplayableUnitaryOperator op)
            {
                if (op?.Data is null)
                {
                    logger.LogError("Asked to encode a displayable unitary operator, but its data was null. This should not happen.");
                    return null;
                }
                return $"Real:\n{op.Data[Ellipsis, 0]}\nImag:\n{op.Data[Ellipsis, 1]}".ToEncodedData();
            }
            else return null;
        }

    }

    public class DisplayableUnitaryOperatorToHtmlEncoder : IResultEncoder
    {
        private readonly ILogger<DisplayableUnitaryOperatorToHtmlEncoder> logger;
        private readonly IConfigurationSource configurationSource;
        public string MimeType => MimeTypes.Html;

        public DisplayableUnitaryOperatorToHtmlEncoder(
            ILogger<DisplayableUnitaryOperatorToHtmlEncoder> logger,
            IConfigurationSource configurationSource
        )
        {
            this.logger = logger;
            this.configurationSource = configurationSource;
        }

        public EncodedData? Encode(object displayable)
        {
            bool IsNearZero(double value) =>
                System.Math.Abs(value) <= configurationSource.TruncationThreshold;

            if (displayable is DisplayableUnitaryOperator op)
            {
                if (op?.Data is null)
                {
                    logger.LogError("Asked to encode a displayable unitary operator, but its data was null. This should not happen.");
                    return null;
                }

                var precision = configurationSource.GetOptionOrDefault(
                    "dump.unitaryPrecision",
                    3
                );

                var outputMatrix = String.Join(
                    " \\\\\n",
                    op.Data.EnumerateOverAxis().Cast<NumSharp.NDArray>().Select(
                        row =>
                            String.Join(" & ",
                                row.EnumerateOverAxis()
                                   .Cast<NumSharp.NDArray>()
                                   .Select(element =>
                                   {
                                       var format = $"{{0:G{precision}}}";
                                       var re = (double)element[0];
                                       var im = (double)element[1];
                                       var reFmt = String.Format(format, re);
                                       var imFmt = String.Format(format, System.Math.Abs(im)) + "i";
                                       if (IsNearZero(re) && IsNearZero(im))
                                       {
                                           return "0";
                                       }
                                       else if (IsNearZero(im))
                                       { 
                                           return reFmt;
                                       }
                                       else if (IsNearZero(re))
                                       {
                                           return im < 0.0 ? $"-{imFmt}" : imFmt;
                                       }
                                       else
                                       {
                                           return $"{reFmt} {(im < 0.0 ? "-" : "+")} {imFmt}";
                                       }
                                   })
                            )
                    )
                );
                return $@"
                    <table>
                        <tr>
                            <th>Qubit IDs</th>
                            <td>{String.Join(", ", op.Qubits.Select(q => q.Id))}
                        </tr>

                        <tr>
                            <th>Unitary representation</th>
                            <td>$$
                                \left(\begin{{matrix}}
                                    {outputMatrix}
                                \end{{matrix}}\right)
                            $$</td>
                        </tr>
                    </table>
                ".ToEncodedData();
            }
            else return null;
        }

    }

}
